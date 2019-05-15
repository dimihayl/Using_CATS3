#include "FemtoBoyzScripts.h"

//#include "math.h"

#include "TH2F.h"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TF1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TColor.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TNtuple.h"
#include "TRandom3.h"

#include "DLM_HistoAnalysis.h"
#include "DLM_SubPads.h"

std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};
std::vector<int> fPalette    = {kWhite, 18, 33, 48, kBlack};
//std::vector<int> fPalette2    = {kWhite, kTeal-9, kGreen-8, kBlack};
std::vector<int> fPalette2    = {kWhite, kRed-1, kRed-9, kBlack};
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyle(bool graypalette=false, bool title=false)
{
      const int NCont = 255;
      gStyle->Reset("Plain");
      gStyle->SetNumberContours(NCont);
      gStyle->SetOptTitle(title);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetOptStat(0);
      if(graypalette) gStyle->SetPalette(8,0);
      else gStyle->SetPalette(1);
      gStyle->SetCanvasColor(10);
      gStyle->SetCanvasBorderMode(0);
      gStyle->SetFrameLineWidth(1);
      gStyle->SetFrameFillColor(kWhite);
      gStyle->SetPadColor(10);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
      gStyle->SetPadBottomMargin(0.15);
      gStyle->SetPadLeftMargin(0.15);
      gStyle->SetHistLineWidth(1);
      gStyle->SetHistLineColor(kRed);
      gStyle->SetFuncWidth(2);
      gStyle->SetFuncColor(kGreen);
      gStyle->SetLineWidth(2);
      gStyle->SetLabelSize(0.045,"xyz");
      gStyle->SetLabelOffset(0.01,"y");
      gStyle->SetLabelOffset(0.01,"x");
      gStyle->SetLabelColor(kBlack,"xyz");
      gStyle->SetTitleSize(0.05,"xyz");
      gStyle->SetTitleOffset(1.25,"y");
      gStyle->SetTitleOffset(1.2,"x");
      gStyle->SetTitleFillColor(kWhite);
      gStyle->SetTextSizePixels(26);
      gStyle->SetTextFont(42);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(kWhite);
      gStyle->SetLegendFont(42);
      gStyle->SetLegendBorderSize(0);

      const int NRGBs = 6;
      Double_t stops[NRGBs];
      for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);

      Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255., 24./255.};
      Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255., 37./255.};
      Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255., 98./255.};
      TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
}
void SetStyle2(bool graypalette=false, bool title=false)
{
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  Int_t palette[5];
  palette[0] = fPalette[0];//kWhite,33,38,kGray; kWhite,18,38,48
  palette[1] = fPalette[1];
  palette[2] = fPalette[2];
  palette[3] = fPalette[3];
  palette[4] = fPalette[4];
  gStyle->SetPalette(5,palette);

//
//  const int NRGBs = 5;
//  Double_t stops[NRGBs];
//  for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);
//  Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255.};
//  Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255.};
//  Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255.};
//  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
}
void SetStyle3(const bool UsePng, const unsigned NumPlotObjects, TAttLine** ObjectList)
{
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(false);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  if(UsePng) gStyle->SetLineWidth(3);
  else gStyle->SetLineWidth(2);
  for(unsigned uObj=0; uObj<NumPlotObjects; uObj++){
    ObjectList[uObj]->SetLineWidth(UsePng?3:2);
  }
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  Int_t palette[5];
  palette[0] = fPalette[0];//kWhite,33,38,kGray
  palette[1] = fPalette[1];
  palette[2] = fPalette[2];
  palette[3] = fPalette[3];
  palette[4] = fPalette[4];
    //TColor** MyColors = new TColor*[4];
    //MyColors[0] = new TColor(5001, 134, 201, 138);
    //MyColors[1] = new TColor(5002, 84, 167, 89);
    //MyColors[2] = new TColor(5003, 0, 147, 154);
    //MyColors[3] = new TColor(5004, 46, 65, 114);
//134,201, 138
//84,167,89
//0 147 154
//46,65,114


  //palette[0] = kWhite;//kWhite,33,38,kGray
  //palette[1] = 5002;
  //palette[2] = 5003;
  //palette[3] = 5004;
  //palette[4] = kBlack;

  gStyle->SetPalette(5,palette);


//
//  const int NRGBs = 5;
//  Double_t stops[NRGBs];
//  for(int i=0; i<NRGBs; ++i) stops[i] = float(i)/(NRGBs-1);
//  Double_t red[NRGBs]   = { 1.,  29./255., 25./255., 27./255., 32./255.};
//  Double_t green[NRGBs] = { 1., 221./255., 160./255., 113./255., 74./255.};
//  Double_t blue[NRGBs] = {  1., 221./255., 184./255., 154./255., 129./255.};
//  TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
}

void SetStyle4(const bool UsePng, const unsigned NumPlotObjects, TAttLine** ObjectList)
{
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(false);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  if(UsePng) gStyle->SetLineWidth(3);
  else gStyle->SetLineWidth(2);
  for(unsigned uObj=0; uObj<NumPlotObjects; uObj++){
    ObjectList[uObj]->SetLineWidth(UsePng?3:2);
  }
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  Int_t palette[4];
  palette[0] = fPalette2[0];//kWhite,33,38,kGray
  palette[1] = fPalette2[1];
  palette[2] = fPalette2[2];
  palette[3] = fPalette2[3];
  gStyle->SetPalette(4,palette);

}

//if Scale, than the result is scaled according to the bin size
TH2F* RebinTH2F(TH2F* hOriginal, const unsigned& RebinX, const unsigned& RebinY, const bool& Scale){
	//hResolutionMatrix->Rebin2D(iRebin);
	//hResolutionMatrix->RebinX(iRebin);
	//hResolutionMatrix->RebinY(iRebin);
	//const unsigned NumBinsX = hResolutionMatrix->GetXaxis()->GetNbinsX();
	//const unsigned NumBinsY = hResolutionMatrix->GetYaxis()->GetNbinsX();
	if(!RebinX||!RebinY) return NULL;
	TH2F* hResult = NULL;
	if(RebinX==1&&RebinY==1) {hResult = (TH2F*)hOriginal->Clone("hResult"); return hResult;}

	unsigned NumBinsX = hOriginal->GetNbinsX();
	unsigned NumBinsY = hOriginal->GetNbinsY();
	unsigned NewBinsX = NumBinsX/unsigned(RebinX);
	unsigned NewBinsY = NumBinsY/unsigned(RebinY);
	double BinWidthX = hOriginal->GetXaxis()->GetBinWidth(1)*double(RebinX);
	double BinWidthY = hOriginal->GetYaxis()->GetBinWidth(1)*double(RebinX);
	double kMinX = hOriginal->GetXaxis()->GetBinLowEdge(1);
	double kMaxX = kMinX + BinWidthX*NewBinsX;
	double kMinY = hOriginal->GetYaxis()->GetBinLowEdge(1);
	double kMaxY = kMinY + BinWidthY*NewBinsY;
//printf("hResolutionMatrix=%p\n",hResolutionMatrix);
//printf(" NumBinsX=%u\n",NumBinsX);
//printf(" NumBinsY=%u\n",NumBinsY);
//printf(" BW=%f\n",hResolutionMatrix->GetXaxis()->GetBinWidth(1));
//printf(" BinWidthX=%f\n",BinWidthX);
//printf(" NewBinsX=%u\n",NewBinsX);
//printf(" unsigned(iRebin*2)=%u\n",unsigned(iRebin*2));
//printf(" kMinX kMaxX=%f %f\n",kMinX,kMaxX);
//printf(" hResolutionMatrix->GetNbins()=%i\n",hResolutionMatrix->GetNbins());
//printf(" hResolutionMatrix(11,11)=%f\n",hResolutionMatrix->GetBinContent(11,11));
	hResult = new TH2F(TString::Format("hResult"),TString::Format("hResult"),NewBinsX,kMinX,kMaxX,NewBinsY,kMinY,kMaxY);
	double BinValue=0;
	double BinError=0;
	double ScaleFactor=1;
	for(unsigned uBinX=0; uBinX<NewBinsX; uBinX++){
		for(unsigned uBinY=0; uBinY<NewBinsX; uBinY++){
			BinValue=0;
			BinError=0;
			for(unsigned SubBinX=0; SubBinX<RebinX; SubBinX++){
				for(unsigned SubBinY=0; SubBinY<RebinY; SubBinY++){
					BinValue += hOriginal->GetBinContent(uBinX*RebinX+SubBinX+1,uBinY*RebinX+SubBinY+1);
					BinError += pow(hOriginal->GetBinError(uBinX*RebinY+SubBinX+1,uBinY*RebinY+SubBinY+1),2.);
				}
			}
			if(Scale){
                ScaleFactor = double(RebinX)*double(RebinY);
			}
			hResult->SetBinContent(uBinX+1,uBinY+1,BinValue/ScaleFactor);
			hResult->SetBinError(uBinX+1,uBinY+1,sqrt(BinError)/ScaleFactor);
		}
	}

    return hResult;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleHisto(TH1 *histo, int marker, int color)
{
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetXaxis()->SetTitleOffset(1.2);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetTitleOffset(1.25);
    histo->SetMarkerSize(1.25);
    if(marker>=0){
        histo->SetLineWidth(2);
        histo->SetMarkerStyle(fMarkers[marker]);
    }
    if(color>=0){
        histo->SetMarkerColor(fColors[color]);
        histo->SetLineColor(fColors[color]);
    }
}

//pL HM
void SetStyleHisto2(TH1 *histo, int marker, int color, double factor=1)
{
    histo->GetXaxis()->SetLabelSize(0.055*factor);
    histo->GetXaxis()->SetTitleSize(0.06*factor);
    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetXaxis()->SetTitleOffset(1.0);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelSize(0.055*factor);
    histo->GetYaxis()->SetTitleSize(0.06*factor);
    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetTitleOffset(1.0/factor);
    histo->SetMarkerSize(1.25);
    if(marker>=0){
        histo->SetLineWidth(2);
        histo->SetMarkerStyle(fMarkers[marker]);
    }
    if(color>=0){
        histo->SetMarkerColor(fColors[color]);
        histo->SetLineColor(fColors[color]);
    }
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *DrawSystematicError(TH1F* histexp,TH1F *histerr,double errorwidth)
{
      //Input are the experimental histogram with statistical errors only and a histogram containing the errors

      const int histbins = histexp->GetNbinsX();

      TGraphErrors *ge_SysError_C2 = new TGraphErrors();
      for(int i=0;i<histbins;i++)
      {
        //if(histexp->GetBinCenter(i+1) > 0.2) continue;
        if(histexp->GetBinCenter(i+1) > 500) continue;
        ge_SysError_C2->SetPoint(i, histexp->GetBinCenter(i+1), histexp->GetBinContent(i+1));
        //printf("ERR at %.2f is %.3f\n",histexp->GetBinCenter(i+1),histerr->GetBinContent(i+1));
        ge_SysError_C2->SetPointError(i, errorwidth, histerr->GetBinContent(i+1));
        //printf();
      }

      return ge_SysError_C2;
}

TGraphErrors *DrawSystematicError_FAST(TH1F* histexp,TH1F *histerr,TF1* ferr,double errorwidth)
{
      //Input are the experimental histogram with statistical errors only and a histogram containing the errors
      const int histbins = histexp->GetNbinsX();
      //TGraph GR_EXP;
      //GR_EXP.SetName("GR_EXP");
      //GR_EXP.Set(histbins);
      //for(int i=0;i<histbins;i++){
      //  GR_EXP.SetPoint(i,histexp->GetBinCenter(i+1),histexp->GetBinContent(i+1));
      //}
    int systbins;
    TGraph GR_SYS;
    GR_SYS.SetName("GR_SYS");
    //take the errors from the histo
    if(histerr){
        systbins = histerr->GetNbinsX();
        GR_SYS.Set(systbins);
        for(int i=0;i<systbins;i++){
            GR_SYS.SetPoint(i,histerr->GetBinCenter(i+1),histerr->GetBinContent(i+1));
        }
    }
    //in the very old output files there was the TF1 fit of the errors, here we read off those
    else if(ferr){
        systbins = histbins;
        GR_SYS.Set(systbins);
        for(int i=0;i<systbins;i++){
            GR_SYS.SetPoint(i,histexp->GetBinCenter(i+1),ferr->Eval(histexp->GetBinCenter(i+1)*0.001));
        }
    }


      TGraphErrors *ge_SysError_C2 = new TGraphErrors();
      for(int i=0;i<histbins;i++)
      {
          double MOM = histexp->GetBinCenter(i+1);
          double VALUE = histexp->GetBinContent(i+1);
          double REL_ERROR = GR_SYS.Eval(MOM);
          double ABS_ERROR = fabs(VALUE*REL_ERROR);
        //if(histexp->GetBinCenter(i+1) > 0.2) continue;
        if(histexp->GetBinCenter(i+1) > 600) continue;
        ge_SysError_C2->SetPoint(i, MOM, VALUE);
        //printf("ERR at %.2f is %.3f\n",histexp->GetBinCenter(i+1),histerr->GetBinContent(i+1));
        ge_SysError_C2->SetPointError(i, errorwidth, ABS_ERROR);
        //printf();
      }

      return ge_SysError_C2;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower, TGraph *grUpper)
{
  TGraphErrors *grFemtoModel = new TGraphErrors();
  grFemtoModel->SetName(grMedian1->GetName());
  double x, yM1, yLo, yUp;
  int count = 0;
  for(int i=0; i< grMedian1->GetN(); ++i) {
    grMedian1->GetPoint(i, x, yM1);
    grLower->GetPoint(i, x, yLo);
    grUpper->GetPoint(i, x, yUp);
    std::vector<float> yAll;
    yAll.push_back(yM1);
    yAll.push_back(yLo);
    yAll.push_back(yUp);
    std::sort(yAll.begin(), yAll.end());
    grFemtoModel->SetPoint(count, x/1000.f, (yAll[2]+yAll[0])/2.f);
    grFemtoModel->SetPointError(count++, 0, (yAll[2]+yAll[0])/2.f - yAll[0]);
  }
  return grFemtoModel;
}

TGraphErrors *FemtoModelFitBandsSimple(TGraph *grLower, TGraph *grUpper)
{
  TGraphErrors *grFemtoModel = new TGraphErrors();
  double x, yM1, yLo, yUp;
  for(int i=0; i< grLower->GetN(); ++i) {
    grLower->GetPoint(i, x, yLo);
    grUpper->GetPoint(i, x, yUp);
    yM1 = (yUp+yLo)*0.5;
    grFemtoModel->SetPoint(i, x, yM1);
//printf("x=%f; yM1=%f\n",x,yM1);
    grFemtoModel->SetPointError(i, 0, (yUp-yLo)*0.5);
  }
  return grFemtoModel;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *convertInGev(TGraphErrors *gr) {
 TGraphErrors *grOut = new TGraphErrors();
 grOut->SetName(gr->GetName());
 double x,y;
 for(int i=0; i< gr->GetN(); ++i) {
   gr->GetPoint(i, x, y);
   grOut->SetPoint(i, x/1000.f, y);
   grOut->SetPointError(i, gr->GetErrorX(i)/1000.f, gr->GetErrorY(i));
 }
 return grOut;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TGraphErrors *convertHistoInGev(TH1F *gr) {
  TGraphErrors *grOut = new TGraphErrors();
  grOut->SetName(gr->GetName());
  for(int i=0; i< gr->GetNbinsX(); ++i) {
    grOut->SetPoint(i, gr->GetBinCenter(i)/1000.f, gr->GetBinContent(i));
    grOut->SetPointError(i, 0, gr->GetBinError(i));
  }
  return grOut;
}

TH1F* Calculate_CF(TH1F* histRE_relK,TH1F* histME_relK, TString CFname,Double_t normleft,Double_t normright, const char* folder, float spinningDepth)
{
  histRE_relK->Sumw2();
  histME_relK->Sumw2();
  TH1F* Hist_CF = (TH1F*)histRE_relK->Clone(CFname.Data());
  if(strcmp(folder, "") == 0) {
    Double_t norm_relK = histRE_relK->Integral(histRE_relK->FindBin(normleft),histRE_relK->FindBin(normright)) / histME_relK->Integral(histME_relK->FindBin(normleft),histME_relK->FindBin(normright));
    Hist_CF->Divide(histRE_relK,histME_relK,1,norm_relK);
  }
  else {
    histME_relK->Scale(1.f/spinningDepth);
    Hist_CF->Divide(histRE_relK,histME_relK);
  }

  return Hist_CF;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F* add_CF(TH1F* hist_CF1, TH1F* hist_CF2,TString HistName)
{
  //Calculate CFs with error weighting
  TH1F* hist_CF_sum = (TH1F*)hist_CF1->Clone(HistName.Data());

  Int_t NBins = hist_CF_sum->GetNbinsX();

  for(Int_t i=0;i<NBins;i++)
  {
    double CF1_val = hist_CF1->GetBinContent(i+1);
    double CF1_err = hist_CF1->GetBinError(i+1);
    double CF2_val = hist_CF2->GetBinContent(i+1);
    double CF2_err = hist_CF2->GetBinError(i+1);

    //average for bin i:
    if(CF1_val != 0. && CF2_val != 0.)
    {
      double CF1_err_weight = 1./TMath::Power(CF1_err,2.);
      double CF2_err_weight = 1./TMath::Power(CF2_err,2.);

      double CF_sum_average = (CF1_err_weight*CF1_val + CF2_err_weight*CF2_val)/(CF1_err_weight+CF2_err_weight);
      double CF_sum_err = 1./TMath::Sqrt(CF1_err_weight+CF2_err_weight);

      hist_CF_sum->SetBinContent(i+1,CF_sum_average);
      hist_CF_sum->SetBinError(i+1,CF_sum_err);
    }
    else if(CF1_val == 0. && CF2_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF2_val);
      hist_CF_sum->SetBinError(i+1,CF2_err);
    }
    else if(CF2_val == 0 && CF1_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF1_val);
      hist_CF_sum->SetBinError(i+1,CF1_err);
    }

  }
  return hist_CF_sum;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SetStyleGraph(TGraph *graph, int marker, int color)
{
  graph->SetMarkerStyle(fMarkers[marker]);
  graph->SetMarkerColor(fColors[color]);
  graph->SetLineColor(fColors[color]);
  graph->GetYaxis()->SetTitleOffset(1.25);
  //graph->SetLineWidth(4);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void plotMorePreliminaries(const TString system)
{

    //const char expfile[] = "~/Results/LHC17p_fast/AnalysisResults.root";
    const TString expfile = system=="pp"?
    "/home/dmihaylov/CernBox/Femto_pp13/data/AnalysisResults_AndiBernieSystME_CkLambdaLambda_0.root":
    "/home/dmihaylov/CernBox/pPb/DataFullest/AnalysisResults_AndiBernieSystME_CkLambdaLambda_0.root";

    const TString MoritaFile = system=="pp"?
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/ComputeMorita_290818/ComputeMoritaPotentials_ALICE_pp_13TeV.root":
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/ComputeMorita_290818/ComputeMoritaPotentials_ALICE_pPb_5TeV.root";

    //int rebin = 1;
    //if(system!="pp") rebin=5;

    gStyle->SetCanvasPreferGL(1);
    const float right = 0.025;
    const float top = 0.025;

    const double energy = system=="pp"?13:5.02; // TeV

    SetStyle();

  const float normleft = 0.2;
  const float normright = 0.4;
  const float spinningDepth = 10.f;
  const char *prefix = "MB";
  const char *addon = "";
  // EXP DATA
  TFile* _file0=TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample6/AnalysisResults.root");
//printf("_file0=%p\n",_file0);
  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon),Results);
  TList* tmpFolder=(TList*)Results->FindObject("Particle0_Particle0");
  tmpFolder=(TList*)Results->FindObject("Particle2_Particle2");
  TH1F* histRE_relK_LL = (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
  TH1F* histME_relK_LL = (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
  if(system=="pp"){
  }
  else{
	//histRE_relK_LL->Rebin(5);
	//histME_relK_LL->Rebin(5);
  }

  tmpFolder=(TList*)Results->FindObject("Particle3_Particle3");
  TH1F* histRE_relK_ALAL = (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
  TH1F* histME_relK_ALAL = (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
  TH1F *hist_CF_LL_ALAL_exp[3];

  hist_CF_LL_ALAL_exp[0] = Calculate_CF(histRE_relK_LL,histME_relK_LL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[1] = Calculate_CF(histRE_relK_ALAL,histME_relK_ALAL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[2] = add_CF(hist_CF_LL_ALAL_exp[0],hist_CF_LL_ALAL_exp[1],"hist_CF_LL_ALAL_exp_sum");

    TString HistoName;
    if(system=="pp") HistoName = "hCkTotNormWeight";
    else HistoName = "hCk_ReweightedMeV_1";

    TFile* fexpfile = new TFile(expfile, "read");
    TH1F* hCkLL_MeV = (TH1F*)fexpfile->Get(HistoName);
    //hCkLL_MeV->Rebin(rebin);
    for(unsigned uBin=0; uBin<20; uBin++){
        hist_CF_LL_ALAL_exp[2]->SetBinContent(uBin+1,hCkLL_MeV->GetBinContent(uBin+1));
        hist_CF_LL_ALAL_exp[2]->SetBinError(uBin+1,hCkLL_MeV->GetBinError(uBin+1));
    }
    hCkLL_MeV = hist_CF_LL_ALAL_exp[2];

    SetStyleHisto(hCkLL_MeV, 0,0);
    hCkLL_MeV->SetMarkerSize(1.5);

    TString input_sys_LL = "/home/dmihaylov/CernBox/Femto_pp13/C2totalsysLL.root";
    TFile* file_sys_LL = new TFile(input_sys_LL,"read");
    TH1F* hist_sys_LL = (TH1F*)file_sys_LL->Get("C2totalsysLL");
    hist_sys_LL->SetLineWidth(2.0);

    TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError(hCkLL_MeV, hist_sys_LL, 0.005);

    TGraph *grFakeSyst = new TGraph();
    grFakeSyst->SetFillColor(fFillColors[0]);
    //grFakeSyst->SetLineColor(kBlack);
    grFakeSyst->SetLineColor(grFakeSyst->GetFillColor());
    grFakeSyst->SetMarkerStyle(20);
    grFakeSyst->SetMarkerSize(1.5);

    TFile* morita = new TFile(MoritaFile,"read");
    auto* grnd46 = convertInGev((TGraphErrors*)morita->Get("gCkTheory_ND46"));
    grnd46->SetFillColor(fColors[2]);
    grnd46->SetLineColor(fColors[2]);
    grnd46->SetLineWidth(3);
    grnd46->SetLineStyle(3);
    //auto* grnsc97f = convertInGev((TGraphErrors*)morita->Get("gCkTheory_NSC97f"));
    //grnsc97f->SetFillColor(fColors[5]);
    //grnsc97f->SetLineColor(fColors[5]);
    //grnsc97f->SetLineWidth(3);
    //grnsc97f->SetLineStyle(2);
    auto* grnf44 = convertInGev((TGraphErrors*)morita->Get("gCkTheory_NF44"));
    grnf44->SetFillColor(fColors[3]);
    grnf44->SetLineColor(fColors[3]);
    grnf44->SetLineWidth(3);
    grnf44->SetLineStyle(2);
    auto* grESC08 = convertInGev((TGraphErrors*)morita->Get("gCkTheory_ESC08"));
    grESC08->SetFillColor(fColors[4]);
    grESC08->SetLineColor(fColors[4]);
    grESC08->SetLineWidth(3);
    auto* grEhime = convertInGev((TGraphErrors*)morita->Get("gCkTheory_Ehime"));
    grEhime->SetFillColor(fColors[2]);
    grEhime->SetLineColor(fColors[2]);
    grEhime->SetLineWidth(3);
    auto* grHKMYY = convertInGev((TGraphErrors*)morita->Get("gCkTheory_HKMYY"));
    grHKMYY->SetFillColor(fColors[5]);
    grHKMYY->SetLineColor(fColors[5]);
    grHKMYY->SetLineWidth(3);
    auto* grNoSI = convertInGev((TGraphErrors*)morita->Get("gCkTheory_NoSI"));
    grNoSI->SetFillColor(kGray+1);
    grNoSI->SetLineColor(kGray+1);
    grNoSI->SetLineWidth(3);
    grNoSI->SetLineStyle(5);

    gStyle->SetErrorX(0.001);

    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.495, 0.875, "ALICE Preliminary");
    if(system=="pp") BeamText.DrawLatex(0.495, 0.825, Form("%s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
    else BeamText.DrawLatex(0.495, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));


    TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
    Can_CF_LL->cd(0);
    Can_CF_LL->SetRightMargin(right);
    Can_CF_LL->SetTopMargin(top);
    Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);
    Tgraph_syserror_LL_ALAL->Draw("ap");
    Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
    Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0.35, 2.);
    Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
    Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.35, 2.);
    grNoSI->Draw("l same");
    grnf44->Draw("l same");
    grEhime->Draw("le3 same");
    grESC08->Draw("le3 same");
    grnd46->Draw("l same");
    grHKMYY->Draw("le3 same");
    //grnsc97f->Draw("l same");
    Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror_LL_ALAL->Draw("2 same");
    hCkLL_MeV->Draw("pe same");
    TLegend *legLL2 = new TLegend(0.5,0.795,0.8,0.86);
    legLL2->SetBorderSize(0);
    legLL2->SetTextFont(42);
    legLL2->SetTextSize(gStyle->GetTextSize()*0.75);
    legLL2->AddEntry(grFakeSyst, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs", "pef");
    legLL2->Draw("same");
    TLegend *legLL = new TLegend(0.5,0.56,0.95,0.76);
    legLL->SetBorderSize(0);
    legLL->SetTextFont(42);
    legLL->SetTextSize(gStyle->GetTextSize()*0.75);
    legLL->SetNColumns(2);
    legLL->AddEntry(grnd46, "ND46", "l");
    //legLL->AddEntry(grnsc97f, "NSC97f", "l");
    legLL->AddEntry(grnf44, "NF44", "l");
    legLL->AddEntry(grEhime, "Ehime", "l");
    legLL->AddEntry(grESC08, "ESC08", "l");
    legLL->AddEntry(grHKMYY, "HKMYY", "l");
    legLL->AddEntry(grNoSI, "Quantum statistics", "l");
    legLL->Draw("same");
    //ref.DrawLatex(0.555, 0.63, "PRC 91 (2015) 024916.");
    if(system=="pp") BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
    else BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
    TDatime DateTime;
    Can_CF_LL->Print(TString::Format("./FemtoBoyzPlots/%i-%i-%i-CF_LL_Models_%s.pdf",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay(),system.Data()));

    delete fexpfile;
    delete file_sys_LL;
}

void plotLambda(const int flag = 0, const TString& Name_Nsigma="hGlobNsigma") {
  SetStyle2();
  //gStyle->SetCanvasPreferGL(1);

  TH1D* hist_emptyhist;
  if(flag==0) hist_emptyhist = new TH1D("hist_emptyhist","; #it{f}_{0}^{#minus1} (fm^{#minus1}); #it{d}_{0} (fm)",50,-2,5);
  if(flag==1) hist_emptyhist = new TH1D("hist_emptyhist","; #it{f}_{0}^{#minus1} (fm^{#minus1}); #it{d}_{0} (fm)",50,-2,5);
  if(flag==2) hist_emptyhist = new TH1D("hist_emptyhist","; #it{f}_{0}^{#minus1} (fm^{#minus1}); #it{d}_{0} (fm)",50,-2,5);
  if(flag==3) hist_emptyhist = new TH1D("hist_emptyhist","; #it{f}_{0}^{#minus1} (fm^{#minus1}); #it{d}_{0} (fm)",50,-0.5,0);

  // Data
  TFile *file2;
  if(flag == 0) file2 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/SRB_111018/Full_7TeV_BL/SystematicsRadiusBinning.root");
  else if(flag == 1) file2 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/Version_270918_2/RT100/CombinedMap.root");
  else if(flag == 2) file2 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/Version_270918_2/RT100/CombinedMap.root");
  else if(flag == 3) file2 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/SRB_111018/Full_7TeV_BL/SystematicsRadiusBinning.root");
  else {printf("The flag is bad!\n"); abort();}

  TH2F *sigma, *ledniSucks,*ledniSucks2;
  bool GetInvLedniSucks = true;
  TFile *file3;
  if (GetInvLedniSucks){
	  if(flag==0) file3 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/SRB_111018/Full_7TeV_BL/Friendship.root");
	  if(flag==1) file3 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/Version_270918_2/RT100/Friendship.root");
	  if(flag==2) file3 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/Version_270918_2/RT100/Friendship.root");
	  if(flag==3) file3 = TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/SRB_111018/Full_7TeV_BL/Friendship.root");
  }

  TH2F* h2Eb;

  // Combined fit
  if(flag == 0) {
    sigma = (TH2F*)file2->Get(Name_Nsigma);
    ledniSucks = (TH2F*)file2->Get("hGlobMinCk");
    if (GetInvLedniSucks) {
      ledniSucks2 = (TH2F*)file3->Get("invLedniSucks");
    }
  }
  // pp 13 TeV
  if(flag == 1) {
    sigma = (TH2F*)file2->Get("hNsigma_1");
    ledniSucks = (TH2F*)file2->Get("hMinCk_1");
  }

  // p#minusPb 5.02 TeV
  if(flag == 2) {
    sigma = (TH2F*)file2->Get("hNsigma_2");
    ledniSucks = (TH2F*)file2->Get("hMinCk_2");
  }

  // Combined fit zoomed in
  if(flag == 3) {
    sigma = (TH2F*)file2->Get(Name_Nsigma);
    ledniSucks = (TH2F*)file2->Get("hGlobMinCk");
    if (GetInvLedniSucks) {
      ledniSucks2 = (TH2F*)file3->Get("invLedniSucks");
    }
    h2Eb = (TH2F*)file2->Get("h2Eb");
    h2Eb->SetLineWidth(2);
    h2Eb->SetLineColor(kBlack);
    h2Eb->SetLineStyle(3);
  }

  SetStyleHisto(sigma);
  sigma->SetMaximum(5);
  sigma->SetMinimum(1);
  sigma->SetTitle(";;;#it{n_{#sigma}}");
  sigma->GetZaxis()->SetLabelFont(hist_emptyhist->GetXaxis()->GetLabelFont());
  sigma->GetZaxis()->SetTitleFont(hist_emptyhist->GetXaxis()->GetTitleFont());
  sigma->GetZaxis()->SetLabelSize(hist_emptyhist->GetXaxis()->GetLabelSize());
  sigma->GetZaxis()->SetTitleSize(hist_emptyhist->GetXaxis()->GetTitleSize());
  TH2F *sigma2 = (TH2F*)sigma->Clone();
  const int nContours = 5;
  double contours[nContours];
    if(Name_Nsigma=="hGlobNsigma"){
      contours[0] = 0;
      contours[1] = 2;
      contours[2] = 3;
      contours[3] = 5;
      contours[4] = 1e6;
    }
    else if(Name_Nsigma=="hGlobChi2Diff"){
      contours[0] = 0;
      contours[1] = 2.3;
      contours[2] = 6.18;
      contours[3] = 11.8;
      contours[4] = 1e6;

      contours[0] = 0;
      contours[1] = sqrt(2.*39.);
      contours[2] = 2.*sqrt(2.*39.);
      contours[3] = 3.*sqrt(2.*39.);
      contours[4] = 1e6;
    }
    else{
      contours[0] = 0;
      contours[1] = 1;
      contours[2] = 3;
      contours[3] = 5;
      contours[4] = 1e6;
    }

  sigma->SetContour(nContours,contours);
  sigma2->SetContour(nContours,contours);

  const int nContoursLedni = 1;
  double contoursLedni[nContoursLedni];
  contoursLedni[0] = 0.00000001;
  ledniSucks->SetContour(nContoursLedni, contoursLedni);
  ledniSucks->SetLineColorAlpha(kGray,0.5);
  ledniSucks->SetLineWidth(2);

  if (GetInvLedniSucks) {
    const int nContoursLedni2 = 5;
    double contoursLedni2[nContoursLedni2];
    contoursLedni2[0] = 0;
    contoursLedni2[1] = 2;
    contoursLedni2[2] = 3;
    contoursLedni2[3] = 5;
    contoursLedni2[4] = 14;
    ledniSucks2->SetFillColor(kBlack);
    //ledniSucks2->SetContour(nContoursLedni2, contoursLedni2);
    //gStyle->SetPalette(5,contoursLedni2);
    //ledniSucks2->SetLineColorAlpha(kGray,0.5);
    //ledniSucks2->SetFillColorAlpha(kGray,0.5);
    ledniSucks2->SetFillStyle(3244);//3663
    ledniSucks2->SetLineColor(kBlack);
  }

  // PRL C02 (2015) 022301.
  TGraphErrors *grStar = new TGraphErrors();
  //https://drupal.star.bnl.gov/STAR/files/starpublications/221/data.html
  grStar->SetPoint(0, -0.91, 8.52);
  grStar->SetPointError(0, 0.31, 2.56);
  SetStyleGraph(grStar, 8, 0);
  grStar->SetMarkerSize(2);
  TGraphAsymmErrors *grStar_syst = new TGraphAsymmErrors();
  grStar_syst->SetPoint(0, -0.91, 8.52);
  grStar_syst->SetPointError(0, 0.56, 0.07, 0.74, 2.09);
  grStar_syst->SetFillColorAlpha(kBlack, 0.4);

//  TGraph *grNagara = new TGraph();
//  grNagara->SetPoint(0, -1.74, 6.45);
//  grNagara->SetPoint(1, -1.3, 6.59);
//  SetStyleGraph(grNagara, 8, 1);
//  grNagara->SetMarkerSize(2);

  // ND
  // M. M. Nagels, T. A. Rijken, and J. J. de Swart, Phys. Rev. D 15, 2547 (1977).
  // Phys. Rev. D 15, 2547 (1977).
  float NDoneOvera0[] = {- 1./4.621, - 1./14.394, 1./10.629, 1./3.483, 1./1.893, 1./1.179, 1./0.764};
  float NDrEff[] = {1.300, 1.633, 2.042, 2.592, 3.389, 4.656, 6.863};
  TGraph *grND = new TGraph(7, NDoneOvera0, NDrEff);
  SetStyleGraph(grND, 2, 2);
  grND->SetLineStyle(3);

  // NF
  // M. M. Nagels, T. A. Rijken, and J. J. de Swart, Phys. Rev. D 20, 1633 (1979).
  // Phys. Rev. D 20, 1633 (1979).
  float NFoneOvera0[] = {- 1./3.659, - 1./23.956, 1./3.960, 1./1.511, 1./0.772, 1./0.406};
  float NFrEff[] = {0.975, 1.258, 1.721, 2.549, 4.271, 8.828};
  TGraph *grNF = new TGraph(6, NFoneOvera0, NFrEff);
  SetStyleGraph(grNF, 3, 3);
  grNF->SetLineStyle(3);

  // NSC89
  // P. M. M. Maessen, T. A. Rijken, and J. J. de Swart, Phys. Rev. C 40, 2226 (1989).
  // Phys. Rev. C 40, 2226 (1989).
  float NSC89oneOvera0[] = {1./0.25, 1./2.1, 1./1.1};
  float NSC89rEff[] = {7.2, 1.9, 3.2};
  TGraph *grNSC89 = new TGraph(3, NSC89oneOvera0, NSC89rEff);
  SetStyleGraph(grNSC89, 4, 4);
  grNSC89->SetLineStyle(2);
  grNSC89->SetMarkerSize(1.4);

  // NSC97
  // T. A. Rijken, V. G. J. Stoks, and Y. Yamamoto, Phys. Rev. C 59, 21 (1999).
  // Phys. Rev. C 59, 21 (1999).
  float NSC97oneOvera0[] = {1./0.329, 1./0.397, 1./0.476, 1./0.401, 1./0.501, 1/0.350};
  float NSC97rEff[] = {12.370, 10.360, 9.130, 1.150, 9.840, 16.330};
  TGraph *grNSC97 = new TGraph(6, NSC97oneOvera0, NSC97rEff);
  SetStyleGraph(grNSC97, 4, 5);
  grNSC97->SetLineStyle(2);
  grNSC97->SetMarkerSize(1.4);

  // Ehime
  // T. Ueda, K. Tominaga, M. Yamaguchi, N. Kijima, D. Okamoto, K. Miyagawa, and T. Yamada, Prog. Theor. Phys. 99, 891 (1998);
  // K. Tominaga, T. Ueda, M. Yamaguchi, N. Kijima, D. Okamoto, K. Miyagawa, and T. Yamada, Nucl. Phys. A 642, 483 (1998).
  // Prog. Theor. Phys. 99, 891 (1998); Nucl. Phys. A 642, 483 (1998).
  TGraph *grEhime = new TGraph();
  grEhime->SetPoint(0, 1./4.21, 2.51);
  SetStyleGraph(grEhime, 0, 2);

  // fss2
  // Y. Fujiwara, Y. Suzuki, and C. Nakamoto, Prog. Part. Nucl. Phys. 58, 439 (2007);
  // Y. Fujiwara, M. Kohno, C. Nakamoto, and Y. Suzuki, Phys. Rev. C 64, 054001 (2001).
  // Part. Nucl. Phys. 58, 439 (2007); Phys. Rev. C 64, 054001 (2001).
  TGraph *grfss2 = new TGraph();
  grfss2->SetPoint(0, 1./0.81, 3.99);
  SetStyleGraph(grfss2, 1, 3);

  // ESC8
  // T. A. Rijken, M. M. Nagels, and Y. Yamamoto, Prog. Theor. Phys. Suppl. 185, 14 (2010).
  // Prog. Theor. Phys. Suppl. 185, 14 (2010).
  TGraph *grESC8 = new TGraph();
  grESC8->SetPoint(0, 1./0.97, 3.86);
  SetStyleGraph(grESC8, 7, 4);
  grESC8->SetMarkerSize(1.4);


  //K. Sasaki and T. Hatsuda (HAL QCD Collaboration), private communication.
  TGraphErrors *grLattice = new TGraphErrors();
  //SetStyleGraph(grLattice, 34, 4);
  grLattice->SetPoint(0, 1.45, 5.16);
  grLattice->SetPointError(0, 0.25, 0.82);
  grLattice->SetMarkerStyle(0);
  grLattice->SetMarkerSize(0);
  grLattice->SetMarkerColor(kRed);
  grLattice->SetLineColor(kRed);
  grLattice->SetLineWidth(2);
  grLattice->SetFillColor(kWhite);

  // HKMYY - NAGARA
  // E. Hiyama, M. Kamimura, T. Motoba, T. Yamada, and Y. Yamamoto, Phys.Rev. C 66, 024007 (2002).
  // Phys.Rev. C 66, 024007 (2002).
  TGraph *grHKMYY = new TGraph();
  grHKMYY->SetPoint(0, 1./0.575, 6.45);
  SetStyleGraph(grHKMYY, 8, 5);
  grHKMYY->SetMarkerSize(2);

  // FG - NAGARA
  // I.N. Filikhin and A. Gal, Nucl. Phys. A 707, 491 (2002).
  // Nucl. Phys. A 707, 491 (2002).
  TGraph *grFG = new TGraph();
  grFG->SetPoint(0, 1./0.77, 6.59);
  SetStyleGraph(grFG, 8, 6);
  grFG->SetMarkerSize(2);

  const float splitFraction = 0.625;
  TCanvas *c = new TCanvas("c", "LambdaLambda", 700./splitFraction, 500);
  TPad *pad1 = new TPad("lambda", "", 0.0,0.,splitFraction,1);
  TPad *pad2 = new TPad("lambda", "", splitFraction,0.,1,1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  pad1->SetTopMargin(0.025);
  pad1->SetRightMargin(0.01);
  if(flag==0){
	hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  }
  if(flag==1){
	hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  }
  if(flag==2){
	hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  }
  if(flag==3){
	hist_emptyhist->GetXaxis()->SetRangeUser(-0.6, 0);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,4.0);
  }

  hist_emptyhist->Draw();
  sigma->Draw("cont3 same");
  sigma->SetLineColor(kGray+3);
  sigma->SetLineWidth(2);
  grStar_syst->Draw("2same");
  grStar->Draw("PZ same");
//  grNagara->Draw("P same");
  grND->Draw("PL same");
  grNF->Draw("PL same");
  grNSC89->Draw("PL same");
  grNSC97->Draw("PL same");
  grEhime->Draw("PL same");
  grfss2->Draw("PL same");
  grESC8->Draw("PL same");
  grHKMYY->Draw("PL same");
  grLattice->Draw("PL same");
  grFG->Draw("PL same");
  //if(flag==3) h2Eb->Draw("cont3,same");

  TGraphErrors *grStarFake = new TGraphErrors();
  grStarFake->SetMarkerStyle(grStar->GetMarkerStyle());
  grStarFake->SetMarkerColor(grStar->GetMarkerColor());
  grStarFake->SetMarkerSize(grStar->GetMarkerSize());
  grStarFake->SetFillColor(fFillColors[0]);
  grStarFake->SetLineColor(fFillColors[0]);

  TGraphErrors *grEntry1Fake = new TGraphErrors();
  grEntry1Fake->SetMarkerColor(fPalette[0]);
  grEntry1Fake->SetMarkerStyle(1);
  grEntry1Fake->SetFillColor(fPalette[0]);
  grEntry1Fake->SetLineColor(kBlack);

  TGraphErrors *grEntry2Fake = new TGraphErrors();
  grEntry2Fake->SetMarkerColor(fPalette[1]);
  grEntry2Fake->SetMarkerStyle(1);
  grEntry2Fake->SetFillColor(fPalette[1]);
  grEntry2Fake->SetLineColor(kBlack);

  TGraphErrors *grEntry3Fake = new TGraphErrors();
  grEntry3Fake->SetMarkerColor(fPalette[2]);
  grEntry3Fake->SetMarkerStyle(1);
  grEntry3Fake->SetFillColor(fPalette[2]);
  grEntry3Fake->SetLineColor(kBlack);

  TGraphErrors *grEntry4Fake = new TGraphErrors();
  grEntry4Fake->SetMarkerColor(fPalette[3]);
  grEntry4Fake->SetMarkerStyle(1);
  grEntry4Fake->SetFillColor(fPalette[3]);
  grEntry4Fake->SetLineColor(kBlack);

  TGraphErrors *grEntry0Fake = new TGraphErrors();
  grEntry0Fake->SetMarkerColor(kWhite);
  grEntry0Fake->SetMarkerStyle(1);
  grEntry0Fake->SetFillColor(kWhite);
  grEntry0Fake->SetLineColor(kWhite);

  pad2->cd();
  pad2->SetTopMargin(0.025);
  TLegend *leg = new TLegend(0.01, 0.18, 0.4, 0.95);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize()*0.85);
  leg->AddEntry(grStarFake, "STAR", "PF");
  leg->AddEntry(grHKMYY, "HKMYY", "P");
  leg->AddEntry(grLattice, "HAL QCD", "L,E0");
  leg->AddEntry(grFG, "FG", "P");
  leg->AddEntry(grND, "ND", "LP");
  leg->AddEntry(grNF, "NF", "LP");
  leg->AddEntry(grNSC89, "NSC89", "LP");
  leg->AddEntry(grNSC97, "NSC97", "LP");
  leg->AddEntry(grEhime, "Ehime", "P");
  leg->AddEntry(grESC8, "ESC08", "P");
  leg->AddEntry(grfss2, "fss2", "P");
  //if(flag==3) leg->AddEntry(h2Eb, "const B_{#Lambda#Lambda}", "L");
  leg->Draw("same");

  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize()*0.68);
  ref.SetNDC(kTRUE);

  const float xAll = 0.35;
  ref.DrawLatex(xAll, 0.895, "PRL C02 (2015) 022301."); // STAR
  ref.DrawLatex(xAll, 0.81828, "Phys.Rev. C 66, 024007 (2002)."); // HKMYY
  ref.DrawLatex(xAll, 0.74162, "Nucl. Phys. A 707, 491 (2002).");  //  FG
  ref.DrawLatex(xAll, 0.66496, "Phys. Rev. D 15, 2547 (1977).");  //  ND
  ref.DrawLatex(xAll, 0.58833, "Phys. Rev. D 20, 1633 (1979).");  //  NF
  ref.DrawLatex(xAll, 0.51164, "Phys. Rev. C 40, 2226 (1989).");  //  NSC89
  ref.DrawLatex(xAll, 0.43599, "Phys. Rev. C 59, 21 (1999).");  //  NSC97
  ref.DrawLatex(xAll, 0.35833, "Nucl. Phys. A 642, 483 (1998).");  //  Ehime Prog. Theor. Phys. 99, 891 (1998);
  ref.DrawLatex(xAll, 0.28166, "Phys. Rev. C 64, 054001 (2001).");  //  fss2 Part. Nucl. Phys. 58, 439 (2007);
  ref.DrawLatex(xAll, 0.205, "Prog. Theor. Phys. Suppl. 185, 14 (2010).");  //  ESC08

  grStar_syst->SetFillColorAlpha(kGray+1, 0.5);
  const float splitFraction2 = 0.8;
  TCanvas *c2 = new TCanvas("c2", "LambdaLambda", 700./splitFraction2, 500);
  TPad *pad12 = new TPad("lambda", "", 0.0,0.,splitFraction2,1);
  TPad *pad22 = new TPad("lambda", "", splitFraction2,0.,1,1);
  pad12->Draw();
  pad22->Draw();
  pad12->cd();
  pad12->SetTopMargin(0.025);
  pad12->SetRightMargin(0.01);
  hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
  hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);

  if(flag==0){
	hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  }
  if(flag==1){
	hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  }
  if(flag==2){
	hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);
  }
  if(flag==3){
	hist_emptyhist->GetXaxis()->SetRangeUser(-0.6, 0);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,4.0);
  }

  hist_emptyhist->Draw();
  sigma2->DrawCopy("col same");
  sigma->Draw("cont3 same");
  sigma->SetLineColor(kWhite);
  sigma->SetLineWidth(1);
  sigma->SetLineStyle(2);

  if (GetInvLedniSucks) {
    ledniSucks2->Draw("BOX same");
  } else {
    ledniSucks->Draw("cont3 same");
  }


  grStar_syst->Draw("2same");
  grStar->Draw("PZ same");
//  grNagara->Draw("P same");
  grND->Draw("PL same");
  grNF->Draw("PL same");
  grNSC89->Draw("PL same");
  grNSC97->Draw("PL same");
  grEhime->Draw("PL same");
  grfss2->Draw("PL same");
  grESC8->Draw("PL same");
  grHKMYY->Draw("PL same");
  grLattice->Draw("PL same");
  grFG->Draw("PL same");
  //if(flag==3) h2Eb->Draw("cont3,same");

  pad12->RedrawAxis();

  TLatex sigmaLabel;
  sigmaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  sigmaLabel.SetTextFont(62);
  sigmaLabel.SetTextColor(kBlack);
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize()*0.8);
  BeamText.SetNDC(kTRUE);

  if(flag == 0) {
//    sigmaLabel.SetTextAngle(32);
//    sigmaLabel.DrawLatex(1.35, 7.4, "1#kern[0.3]{#sigma}");
    if(Name_Nsigma=="hGlobNsigma"){
    sigmaLabel.SetTextAngle(-50);
    sigmaLabel.DrawLatex(0.50, 1.4, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-55);
    sigmaLabel.DrawLatex(0.20, 1.2, "5#kern[0.2]{#sigma}");

    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.6, 1.95, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.6, 3.3, "5#kern[0.2]{#sigma}");
    }


    BeamText.SetTextSize(gStyle->GetTextSize()*0.68);
    BeamText.DrawLatex(0.73, 0.36, "ALICE");
    BeamText.DrawLatex(0.73, 0.315, "pp #sqrt{#it{s}} = 7 TeV");
    //BeamText.DrawLatex(0.73, 0.315, "ALICE");
    BeamText.DrawLatex(0.73, 0.27, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.73, 0.225, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    BeamText.DrawLatex(0.73, 0.18, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs");
  }
  if(flag == 1) {
    sigmaLabel.SetTextAngle(12.5);
    sigmaLabel.DrawLatex(1.35, 7.5, "1#kern[0.3]{#sigma}");
    sigmaLabel.SetTextAngle(-50);
    sigmaLabel.DrawLatex(0.7, 1.5, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-65);
    sigmaLabel.DrawLatex(0.35, 1.5, "5#kern[0.2]{#sigma}");

    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.2, 1.5, "3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.05, 2.7, "5#kern[0.2]{#sigma}");

    BeamText.DrawLatex(0.6, 0.305, "ALICE");
    BeamText.DrawLatex(0.6, 0.255, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.6, 0.205, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs");
  }
  if(flag == 2) {
    sigmaLabel.SetTextAngle(320);
    sigmaLabel.DrawLatex(-1.74,1.16,"3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-35);
    sigmaLabel.DrawLatex(-1.74,3.41,"5#kern[0.2]{#sigma}");

    sigmaLabel.SetTextAngle(300);
    sigmaLabel.DrawLatex(0.41,1.25,"3#kern[0.2]{#sigma}");
    sigmaLabel.SetTextAngle(-65);
    sigmaLabel.DrawLatex(0.15,1.25, "5#kern[0.2]{#sigma}");

    BeamText.DrawLatex(0.6, 0.305, "ALICE");
    BeamText.SetTextSize(gStyle->GetTextSize()*0.65);
    BeamText.DrawLatex(0.6, 0.255, "p#minusPb #sqrt{#it{s_{NN}}} = 5.02 TeV");
    BeamText.DrawLatex(0.6, 0.205, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs");
  }
  if(flag == 3) {
//    sigmaLabel.SetTextAngle(32);
//    sigmaLabel.DrawLatex(1.35, 7.4, "1#kern[0.3]{#sigma}");

    sigmaLabel.SetTextAngle(60);
    sigmaLabel.SetTextSize(gStyle->GetTextSize()*0.7);
    sigmaLabel.SetTextFont(42);

    //sigmaLabel.DrawLatex(-0.210, 0.12, "B_{#Lambda#Lambda}=1.9 MeV");
    //sigmaLabel.DrawLatex(-0.280, 0.12, "B_{#Lambda#Lambda}=2.5 MeV");
    //sigmaLabel.DrawLatex(-0.335, 0.12, "B_{#Lambda#Lambda}=3.8 MeV");

    BeamText.SetTextSize(gStyle->GetTextSize()*0.68);
    BeamText.DrawLatex(0.18, 0.41, "ALICE");
    BeamText.DrawLatex(0.18, 0.365, "pp #sqrt{#it{s}} = 7 TeV");
    //BeamText.DrawLatex(0.18, 0.365, "ALICE");
    BeamText.DrawLatex(0.18, 0.32, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.18, 0.275, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    BeamText.DrawLatex(0.18, 0.23, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs");
  }

  pad22->cd();
  pad22->SetTopMargin(0.025);
  TLegend *leg2 ;
  if(flag!=3) leg2 = new TLegend(0.01, 0.01, 0.99, 0.97);
  else leg2 = new TLegend(0.01, 0.31, 0.99, 0.97);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(gStyle->GetTextSize()*2.25);
//  leg2->AddEntry(grEntry3Fake, "Exclusion", "");
//  leg2->AddEntry(grEntry3Fake, "#sigma < 3", "PF");
  if(flag!=3){
    if(Name_Nsigma=="hGlobNsigma"){
        leg2->AddEntry(grEntry1Fake, "#sigma < 2", "PF");
        leg2->AddEntry(grEntry2Fake, "2 < #sigma < 3", "PF");
        leg2->AddEntry(grEntry3Fake, "3 < #sigma < 5", "PF");
        leg2->AddEntry(grEntry4Fake, "#sigma >5 ", "PF");
    }
    else if(Name_Nsigma=="hGlobChi2Diff"){
        leg2->AddEntry(grEntry1Fake, TString::Format("#Delta#chi^{2}<%.1f",contours[1]), "PF");
        leg2->AddEntry(grEntry2Fake, TString::Format("%.1f<#Delta#chi^{2}<%.1f",contours[1],contours[2]), "PF");
        leg2->AddEntry(grEntry3Fake, TString::Format("%.1f<#Delta#chi^{2}<%.1f",contours[2],contours[3]), "PF");
        leg2->AddEntry(grEntry4Fake, TString::Format("#Delta#chi^{2}>%.1f",contours[3]), "PF");
    }
    else{
        leg2->AddEntry(grEntry1Fake, "#sigma < 1", "PF");
        leg2->AddEntry(grEntry2Fake, "1 < #sigma < 3", "PF");
        leg2->AddEntry(grEntry3Fake, "3 < #sigma < 5", "PF");
        leg2->AddEntry(grEntry4Fake, "#sigma >5 ", "PF");
    }

	leg2->AddEntry(ledniSucks2, "Unphys. C(k*)", "PF");
	leg2->AddEntry(grEntry0Fake, " ", "");
	leg2->AddEntry(grStarFake, "STAR", "PF");
	leg2->AddEntry(grLattice, "HAL QCD", "L,E0");
	leg2->AddEntry(grHKMYY, "HKMYY", "P");
	leg2->AddEntry(grFG, "FG", "P");
	leg2->AddEntry(grND, "ND", "LP");
	leg2->AddEntry(grNF, "NF", "LP");
	leg2->AddEntry(grNSC89, "NSC89", "LP");
	leg2->AddEntry(grNSC97, "NSC97", "LP");
	leg2->AddEntry(grEhime, "Ehime", "P");
	leg2->AddEntry(grESC8, "ESC08", "P");
	leg2->AddEntry(grfss2, "fss2", "P");
  }
  else{
    if(Name_Nsigma=="hGlobNsigma"){
        leg2->AddEntry(grEntry1Fake, "#sigma < 2", "PF");
        leg2->AddEntry(grEntry2Fake, "2 < #sigma < 3", "PF");
        leg2->AddEntry(grEntry3Fake, "3 < #sigma < 5", "PF");
        leg2->AddEntry(grEntry4Fake, "#sigma >5 ", "PF");
    }
    else if(Name_Nsigma=="hGlobChi2Diff"){
        leg2->AddEntry(grEntry1Fake, TString::Format("#Delta#chi^{2}<%.1f",contours[1]), "PF");
        leg2->AddEntry(grEntry2Fake, TString::Format("%.1f<#Delta#chi^{2}<%.1f",contours[1],contours[2]), "PF");
        leg2->AddEntry(grEntry3Fake, TString::Format("%.1f<#Delta#chi^{2}<%.1f",contours[2],contours[3]), "PF");
        leg2->AddEntry(grEntry4Fake, TString::Format("#Delta#chi^{2}>%.1f",contours[3]), "PF");
    }
    else{
        leg2->AddEntry(grEntry1Fake, "#sigma < 1", "PF");
        leg2->AddEntry(grEntry2Fake, "1 < #sigma < 3", "PF");
        leg2->AddEntry(grEntry3Fake, "3 < #sigma < 5", "PF");
        leg2->AddEntry(grEntry4Fake, "#sigma >5 ", "PF");
    }
	leg2->AddEntry(ledniSucks2, "Unphys. C(k*)", "PF");
	leg2->AddEntry(grEntry0Fake, " ", "");
	leg2->AddEntry(grND, "ND 46,48", "LP");
	leg2->AddEntry(grNF, "NF 42,44", "LP");
	//leg2->AddEntry(h2Eb, "const B_{#Lambda#Lambda}", "L");
  }

  leg2->Draw("same");

  TDatime DateTime;
  if(flag == 0) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_global.pdf",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));
  if(flag == 1) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_13TeV.pdf",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));
  if(flag == 2) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_pPb502TeV.pdf",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));
  if(flag == 3) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_globalzoom.pdf",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));

c2->cd(0); c2->SetCanvasSize(1920, 1080); //c2->SetMargin(0.13,0.05,0.2,0.05);//lrbt
  if(flag == 0) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_global.png",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));
  if(flag == 1) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_13TeV.png",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));
  if(flag == 2) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_pPb502TeV.png",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));
  if(flag == 3) c2->SaveAs(TString::Format("./FemtoBoyzPlots/%i-%i-%i-LambdaLambda_globalzoom.png",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay()));



}


//we plot the exclusion plot 1/f0 vs d0, not zoomed
void PlotDimiExclusion_ver1(const TString FolderName, const TString HistoNameEP, const TString FileNameLED, const TString HistoNameLED,
                            const unsigned RebinX, const unsigned RebinY,
                            const TString DataSets){

    gStyle->SetCanvasPreferGL(1);

    const TString OutputFileName = "PlotDimiExclusion_ver1.root";

    //histo with the data
    const TString InputFileName = "FinalEP.root";
    //const TString HistoNameEP = "hFinalEP_Confidence";
    TFile* FileFinalEP = new TFile(FolderName+InputFileName,"read");
    TH2F* hFinalEP = (TH2F*)FileFinalEP->Get(HistoNameEP);
    TH2F* hFinalEP_Rebinned = RebinTH2F(hFinalEP,RebinX,RebinY,true);

    TFile* FileLED = new TFile(FileNameLED,"read");
    TH2F* hNegativeLednicky = (TH2F*)FileLED->Get(HistoNameLED);
    TH2F* hNegativeLednicky_Rebinned = RebinTH2F(hNegativeLednicky,1,1,true);
    for(unsigned uBinX=0; uBinX<hNegativeLednicky_Rebinned->GetNbinsX(); uBinX++){
        for(unsigned uBinY=0; uBinY<hNegativeLednicky_Rebinned->GetNbinsY(); uBinY++){
            hNegativeLednicky_Rebinned->SetBinContent(uBinX,uBinY,
                            hNegativeLednicky_Rebinned->GetBinContent(uBinX,uBinY)?1:0);
        }
    }

    hNegativeLednicky_Rebinned->SetFillColor(kGray+3);
    hNegativeLednicky_Rebinned->SetFillStyle(3344);//3663,3244
    hNegativeLednicky_Rebinned->SetLineColor(kGray+3);
    hNegativeLednicky_Rebinned->SetLineWidth(1.5);
    hNegativeLednicky_Rebinned->SetLineStyle(1);

    const unsigned NumContours = 5;
    double* nSigma_Contours = new double [NumContours];
    nSigma_Contours[0] = 0;
    nSigma_Contours[1] = 1;
    nSigma_Contours[2] = 2;
    nSigma_Contours[3] = 3;
    nSigma_Contours[4] = 100;
    hFinalEP_Rebinned->SetContour(NumContours,nSigma_Contours);
//DO THAT!
    TH2F* hLedniFails;

    const unsigned NumModels=11;
    TString* ModelName = new TString [NumModels];
    unsigned Counter=0;
    ModelName[Counter++] = "STAR";//0
    ModelName[Counter++] = "HAL QCD";
    ModelName[Counter++] = "HKMYY";//2
    ModelName[Counter++] = "FG";
    ModelName[Counter++] = "ND";//4
    ModelName[Counter++] = "NF";
    ModelName[Counter++] = "NSC89";//6
    ModelName[Counter++] = "NSC97";
    ModelName[Counter++] = "Ehime";//8
    ModelName[Counter++] = "ESC08";
    ModelName[Counter++] = "fss2";//10
    if(Counter!=NumModels){
        printf("\033[1;33mWARNING:\033[0m Counter!=NumModels (%u!=%u)\n",Counter,NumModels);
    }

    unsigned* NumModelPoints = new unsigned [NumModels];
    Counter=0;
    NumModelPoints[Counter++]=1;//0
    NumModelPoints[Counter++]=1;
    NumModelPoints[Counter++]=1;//2
    NumModelPoints[Counter++]=1;
    NumModelPoints[Counter++]=7;//4
    NumModelPoints[Counter++]=6;
    NumModelPoints[Counter++]=3;//6
    NumModelPoints[Counter++]=6;
    NumModelPoints[Counter++]=1;//8
    NumModelPoints[Counter++]=1;
    NumModelPoints[Counter++]=1;//10

    TGraphErrors* grModel = new TGraphErrors[NumModels];

    for(unsigned uMod=0; uMod<NumModels; uMod++){
        //grModel[uMod].SetLineWidth(9);
        grModel[uMod].SetMarkerSize(1);
    }

    //STAR
    //https://drupal.star.bnl.gov/STAR/files/starpublications/221/data.html
    grModel[0].Set(1);
    grModel[0].SetName("grSTAR");
    grModel[0].SetPoint(0, -0.91, 8.52);
    grModel[0].SetPointError(0, 0.31, 2.56);
    SetStyleGraph(&grModel[0], 8, 0);
    grModel[0].SetMarkerStyle(29);
    //grModel[0].SetMarkerColor(kBlue+2);
    grModel[0].SetMarkerSize(2.9*1);
    //grModel[0].SetLineColor(kBlue+2);
    grModel[0].SetFillColorAlpha(kGray+2, 0.6);
    TGraphAsymmErrors *grStar_syst = new TGraphAsymmErrors();
    grStar_syst->SetPoint(0, -0.91, 8.52);
    grStar_syst->SetPointError(0, 0.56, 0.07, 0.74, 2.09);
    //grStar_syst->SetFillColorAlpha(kBlack, 0.4);
    grStar_syst->SetFillColorAlpha(kGray+2, 0.6);
    //grStar_syst->SetFillColor(kGray+1);
    //grStar_syst->SetLineColor(kBlack);


    //HAL QCD
    //K. Sasaki and T. Hatsuda (HAL QCD Collaboration), private communication.
    //SetStyleGraph(&grModel[1], 34, 4);
    grModel[1].Set(1);
    grModel[1].SetName("grHALQCD");
    grModel[1].SetPoint(0, 1.45, 5.16);
    grModel[1].SetPointError(0, 0.25, 0.82);
    grModel[1].SetMarkerStyle(0);
    grModel[1].SetMarkerSize(0);
    grModel[1].SetMarkerColor(kRed);
    grModel[1].SetLineColor(kRed);
    //grModel[1].SetLineWidth(3);
    grModel[1].SetFillColor(kWhite);

    // HKMYY - NAGARA
    // E. Hiyama, M. Kamimura, T. Motoba, T. Yamada, and Y. Yamamoto, Phys.Rev. C 66, 024007 (2002).
    // Phys.Rev. C 66, 024007 (2002).
    grModel[2].Set(1);
    grModel[2].SetName("grHKMYY");
    grModel[2].SetPoint(0, 1./0.575, 6.45);
    SetStyleGraph(&grModel[2], 8, 5);
    grModel[2].SetMarkerSize(2.6*1);

    // FG - NAGARA
    // I.N. Filikhin and A. Gal, Nucl. Phys. A 707, 491 (2002).
    // Nucl. Phys. A 707, 491 (2002).
    grModel[3].Set(1);
    grModel[3].SetName("grFG");
    grModel[3].SetPoint(0, 1./0.77, 6.59);
    SetStyleGraph(&grModel[3], 8, 6);
    grModel[3].SetMarkerSize(2.6*1);

    // ND
    // M. M. Nagels, T. A. Rijken, and J. J. de Swart, Phys. Rev. D 15, 2547 (1977).
    // Phys. Rev. D 15, 2547 (1977).
    float NDoneOvera0[7] = {-1./4.621, -1./14.394, 1./10.629, 1./3.483, 1./1.893, 1./1.179, 1./0.764};
    float NDrEff[7] = {1.300, 1.633, 2.042, 2.592, 3.389, 4.656, 6.863};
    grModel[4].Set(7);
    grModel[4].SetName("grND");
    for(unsigned uTmp=0; uTmp<7; uTmp++) grModel[4].SetPoint(uTmp,NDoneOvera0[uTmp],NDrEff[uTmp]);
    SetStyleGraph(&grModel[4], 2, 2);
    grModel[4].SetLineStyle(3);

    // NF
    // M. M. Nagels, T. A. Rijken, and J. J. de Swart, Phys. Rev. D 20, 1633 (1979).
    // Phys. Rev. D 20, 1633 (1979).
    float NFoneOvera0[6] = {-1./3.659, -1./23.956, 1./3.960, 1./1.511, 1./0.772, 1./0.406};
    float NFrEff[6] = {0.975, 1.258, 1.721, 2.549, 4.271, 8.828};
    grModel[5].Set(6);
    grModel[5].SetName("grNF");
    for(unsigned uTmp=0; uTmp<6; uTmp++) grModel[5].SetPoint(uTmp,NFoneOvera0[uTmp],NFrEff[uTmp]);
    SetStyleGraph(&grModel[5],3,3);
    grModel[5].SetLineStyle(3);

    // NSC89
    // P. M. M. Maessen, T. A. Rijken, and J. J. de Swart, Phys. Rev. C 40, 2226 (1989).
    // Phys. Rev. C 40, 2226 (1989).
    float NSC89oneOvera0[3] = {1./0.25, 1./2.1, 1./1.1};
    float NSC89rEff[3] = {7.2, 1.9, 3.2};
    grModel[6].Set(3);
    grModel[6].SetName("grNSC89");
    for(unsigned uTmp=0; uTmp<3; uTmp++) grModel[6].SetPoint(uTmp,NSC89oneOvera0[uTmp],NSC89rEff[uTmp]);
    SetStyleGraph(&grModel[6],4,4);
    grModel[6].SetLineStyle(2);
    grModel[6].SetMarkerSize(1.4*1);

    // NSC97
    // T. A. Rijken, V. G. J. Stoks, and Y. Yamamoto, Phys. Rev. C 59, 21 (1999).
    // Phys. Rev. C 59, 21 (1999).
    float NSC97oneOvera0[6] = {1./0.329, 1./0.397, 1./0.476, 1./0.401, 1./0.501, 1/0.350};
    float NSC97rEff[6] = {12.370, 10.360, 9.130, 1.150, 9.840, 16.330};
    grModel[7].Set(6);
    grModel[7].SetName("grNSC97");
    for(unsigned uTmp=0; uTmp<6; uTmp++) grModel[7].SetPoint(uTmp,NSC97oneOvera0[uTmp],NSC97rEff[uTmp]);
    SetStyleGraph(&grModel[7], 4, 5);
    grModel[7].SetLineStyle(2);
    grModel[7].SetMarkerSize(1.4*1);

    // Ehime
    // T. Ueda, K. Tominaga, M. Yamaguchi, N. Kijima, D. Okamoto, K. Miyagawa, and T. Yamada, Prog. Theor. Phys. 99, 891 (1998);
    // K. Tominaga, T. Ueda, M. Yamaguchi, N. Kijima, D. Okamoto, K. Miyagawa, and T. Yamada, Nucl. Phys. A 642, 483 (1998).
    // Prog. Theor. Phys. 99, 891 (1998); Nucl. Phys. A 642, 483 (1998).
    grModel[8].Set(1);
    grModel[8].SetName("grEhime");
    grModel[8].SetPoint(0, 1./4.21, 2.51);
    SetStyleGraph(&grModel[8], 0, 2);

    // ESC8
    // T. A. Rijken, M. M. Nagels, and Y. Yamamoto, Prog. Theor. Phys. Suppl. 185, 14 (2010).
    // Prog. Theor. Phys. Suppl. 185, 14 (2010).
    grModel[9].Set(1);
    grModel[9].SetName("grESC8");
    grModel[9].SetPoint(0, 1./0.97, 3.86);
    SetStyleGraph(&grModel[9], 7, 4);
    grModel[9].SetMarkerSize(1.4*1);

    // fss2
    // Y. Fujiwara, Y. Suzuki, and C. Nakamoto, Prog. Part. Nucl. Phys. 58, 439 (2007);
    // Y. Fujiwara, M. Kohno, C. Nakamoto, and Y. Suzuki, Phys. Rev. C 64, 054001 (2001).
    // Part. Nucl. Phys. 58, 439 (2007); Phys. Rev. C 64, 054001 (2001).
    grModel[10].Set(1);
    grModel[10].SetName("grfss2");
    grModel[10].SetPoint(0, 1./0.81, 3.99);
    SetStyleGraph(&grModel[10], 1, 3);



    TGraphErrors *grEntry1Fake = new TGraphErrors();
    grEntry1Fake->SetMarkerColor(fPalette[0]);
    grEntry1Fake->SetMarkerStyle(1);
    grEntry1Fake->SetFillColor(fPalette[0]);
    grEntry1Fake->SetLineColor(kBlack);

    TGraphErrors *grEntry2Fake = new TGraphErrors();
    grEntry2Fake->SetMarkerColor(fPalette[1]);
    grEntry2Fake->SetMarkerStyle(1);
    grEntry2Fake->SetFillColor(fPalette[1]);
    grEntry2Fake->SetLineColor(kBlack);

    TGraphErrors *grEntry3Fake = new TGraphErrors();
    grEntry3Fake->SetMarkerColor(fPalette[2]);
    grEntry3Fake->SetMarkerStyle(1);
    grEntry3Fake->SetFillColor(fPalette[2]);
    grEntry3Fake->SetLineColor(kBlack);

    TGraphErrors *grEntry4Fake = new TGraphErrors();
    grEntry4Fake->SetMarkerColor(fPalette[3]);
    grEntry4Fake->SetMarkerStyle(1);
    grEntry4Fake->SetFillColor(fPalette[3]);
    grEntry4Fake->SetLineColor(kBlack);

    TFile* OutputFile = new TFile(FolderName+OutputFileName,"recreate");

    SetStyle3(true,0,0);
    const float splitFraction = 0.8;
    TCanvas *cLamLamExclusion = new TCanvas("cLamLamExclusion", "LambdaLambda", 960, 540);
    TPad *pad1 = new TPad("lambda", "", 0.0,0.,splitFraction,1);
    TPad *pad2 = new TPad("lambda", "", splitFraction,0.,1,1);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    pad1->SetTopMargin(0.025);
    pad1->SetRightMargin(0.01);

    pad2->cd();
    //cLamLamExclusion->cd(0);
    TLegend* leg2 = new TLegend(0.01, 0.71, 0.99, 0.97);//lbrt
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(gStyle->GetTextSize()*2.25);
    leg2->AddEntry(grEntry1Fake, TString::Format("#it{n#sigma} < %.0f",nSigma_Contours[1]), "PF");
    leg2->AddEntry(grEntry2Fake, TString::Format("%.0f < #it{n#sigma} < %.0f",nSigma_Contours[1],nSigma_Contours[2]), "PF");
    leg2->AddEntry(grEntry3Fake, TString::Format("%.0f < #it{n#sigma} < %.0f",nSigma_Contours[2],nSigma_Contours[3]), "PF");
    leg2->AddEntry(grEntry4Fake, TString::Format("#it{n#sigma} > %.0f",nSigma_Contours[3]), "PF");
    leg2->AddEntry(hNegativeLednicky_Rebinned, "Unphys. #it{C(k*)}", "PF");
    leg2->Draw("same");

    TLegend* leg3 = new TLegend(0.01, 0.06, 0.99, 0.66);//lbrt
    leg3->SetBorderSize(0);
    leg3->SetTextFont(42);
    leg3->SetTextSize(gStyle->GetTextSize()*2.25);
    for(unsigned uMod=0; uMod<NumModels; uMod++){
        if(ModelName[uMod]=="HAL QCD") leg3->AddEntry(&grModel[uMod], TString::Format("%s",ModelName[uMod].Data()), "L,E0");
        else if(ModelName[uMod]=="STAR") leg3->AddEntry(&grModel[uMod], TString::Format("%s",ModelName[uMod].Data()));
        else leg3->AddEntry(&grModel[uMod], TString::Format("%s",ModelName[uMod].Data()), (NumModelPoints[uMod]>1)?TString("LP"):TString("P"));
    }
    leg3->Draw("same");

    pad1->cd();

    TH1D* hist_emptyhist;
    //hist_emptyhist = new TH1D("hist_emptyhist","; #it{f}_{0}^{#minus1} (fm^{#minus1}); #it{d}_{0} (fm)",50,-2,5);
    hist_emptyhist = new TH1D("hist_emptyhist","; (fm^{#minus1}); #it{d}_{0} (fm)",50,-2,5);

    hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);

	hist_emptyhist->Draw();
	hFinalEP_Rebinned->Draw("same,col,hist");
	hist_emptyhist->Draw("same,axis");
	grStar_syst->Draw("2same");
	//grStar_syst->Draw("BOX,same");
	hNegativeLednicky_Rebinned->Draw("BOX,same");
	for(unsigned uMod=NumModels-1; uMod>=0; uMod--){
        grModel[uMod].Draw("PL,same");
        if(uMod==0) break;
	}
	//for(unsigned uMod=0; uMod<NumModels; uMod++){
    //    grModel[uMod].Draw("PL,same");
	//}

	//leg2->Draw("same");

    TLatex ScatLenText;
    ScatLenText.SetTextSize(hist_emptyhist->GetXaxis()->GetTitleSize());
    ScatLenText.SetNDC(kTRUE);
    ScatLenText.SetTextFont(12);
    ScatLenText.DrawLatex(0.875,0.042,"f_{#it{0}}^{#it{#minus1}}");
    ScatLenText.Draw();

    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.8);
    BeamText.SetNDC(kTRUE);

    BeamText.SetTextSize(gStyle->GetTextSize()*0.68);
    BeamText.DrawLatex(0.75, 0.36, "ALICE");
    BeamText.DrawLatex(0.75, 0.315, "pp #sqrt{#it{s}} = 7 TeV");
    //BeamText.DrawLatex(0.73, 0.315, "ALICE");
    BeamText.DrawLatex(0.75, 0.27, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.75, 0.225, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    BeamText.DrawLatex(0.75, 0.18, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs");
    BeamText.Draw();

	//pad1->RedrawAxis();
	//pad2->RedrawAxis();
	cLamLamExclusion->Write();

	const unsigned NumPlotObj = NumModels;
	TAttLine** ListOfObjects = new TAttLine* [NumPlotObj];
	for(unsigned uMod=0; uMod<NumModels; uMod++){
        ListOfObjects[uMod] = &grModel[uMod];
	}

	SetStyle3(false,NumPlotObj,ListOfObjects);

	cLamLamExclusion->SaveAs(FolderName+"cLamLamExclusion.pdf");
	SetStyle3(true,NumPlotObj,ListOfObjects);
	cLamLamExclusion->SaveAs(FolderName+"cLamLamExclusion.png");

    delete [] ModelName;
    delete [] NumModelPoints;
    delete [] nSigma_Contours;
    delete [] grModel;
    delete [] ListOfObjects;
    delete grStar_syst;
    delete hist_emptyhist;
    delete pad1;
    delete pad2;
    delete leg2;
    delete leg3;
    delete grEntry1Fake;
    delete grEntry2Fake;
    delete grEntry3Fake;
    delete grEntry4Fake;
    delete hFinalEP_Rebinned;
    delete hNegativeLednicky_Rebinned;
    delete cLamLamExclusion;
    delete FileFinalEP;
    delete FileLED;
    delete OutputFile;
}

//we plot the exclusion plot 1/f0 vs d0, not zoomed
void PlotDimiChi2Exclusion_ver1(const TString FolderName, const TString InputFileName, const TString HistoNameEP,
                            const unsigned RebinX, const unsigned RebinY,
                            const TString DataSets){

    const TString OutputFileName = "PlotDimiChi2Exclusion_ver1.root";

    //const TString HistoNameEP = "hFinalEP_Confidence";
    TFile* FileFinalEP = new TFile(FolderName+InputFileName,"read");
    TH2F* hFinalEP = (TH2F*)FileFinalEP->Get(HistoNameEP);
    TH2F* hFinalEP_Rebinned = RebinTH2F(hFinalEP,RebinX,RebinY,true);


    const unsigned NumContours = 5;
    double* nSigma_Contours = new double [NumContours];
    nSigma_Contours[0] = 0;
    nSigma_Contours[1] = 0.1;
    nSigma_Contours[2] = 0.5;
    nSigma_Contours[3] = 2.5;
    nSigma_Contours[4] = 100;
    //hFinalEP_Rebinned->SetContour(NumContours,nSigma_Contours);

    double MinimumValue;
    double TempValue;
    MinimumValue = hFinalEP_Rebinned->GetMinimum();
    const double Ndof = 61;
    for(unsigned uBinX=1; uBinX<=hFinalEP_Rebinned->GetNbinsX(); uBinX++){
        for(unsigned uBinY=1; uBinY<=hFinalEP_Rebinned->GetNbinsY(); uBinY++){
            TempValue = hFinalEP_Rebinned->GetBinContent(uBinX,uBinY);
            TempValue += MinimumValue;
            //if(TempValue>nSigma_Contours[3]) TempValue = 0.5*(nSigma_Contours[4]+nSigma_Contours[3]);
            hFinalEP_Rebinned->SetBinContent(uBinX,uBinY,TempValue*Ndof);
        }
    }


    TGraphErrors *grEntry1Fake = new TGraphErrors();
    grEntry1Fake->SetMarkerColor(fPalette[0]);
    grEntry1Fake->SetMarkerStyle(1);
    grEntry1Fake->SetFillColor(fPalette[0]);
    grEntry1Fake->SetLineColor(kBlack);

    TGraphErrors *grEntry2Fake = new TGraphErrors();
    grEntry2Fake->SetMarkerColor(fPalette[1]);
    grEntry2Fake->SetMarkerStyle(1);
    grEntry2Fake->SetFillColor(fPalette[1]);
    grEntry2Fake->SetLineColor(kBlack);

    TGraphErrors *grEntry3Fake = new TGraphErrors();
    grEntry3Fake->SetMarkerColor(fPalette[2]);
    grEntry3Fake->SetMarkerStyle(1);
    grEntry3Fake->SetFillColor(fPalette[2]);
    grEntry3Fake->SetLineColor(kBlack);

    TGraphErrors *grEntry4Fake = new TGraphErrors();
    grEntry4Fake->SetMarkerColor(fPalette[3]);
    grEntry4Fake->SetMarkerStyle(1);
    grEntry4Fake->SetFillColor(fPalette[3]);
    grEntry4Fake->SetLineColor(kBlack);

    TFile* OutputFile = new TFile(FolderName+OutputFileName,"recreate");

    SetStyle3(true,0,0);
    const float splitFraction = 0.8;
    TCanvas *cLamLamExclusion = new TCanvas("cLamLamExclusion", "LambdaLambda", 960, 540);
    cLamLamExclusion->cd(0);

    TH1D* hist_emptyhist;
    //hist_emptyhist = new TH1D("hist_emptyhist","; #it{f}_{0}^{#minus1} (fm^{#minus1}); #it{d}_{0} (fm)",50,-2,5);
    hist_emptyhist = new TH1D("hist_emptyhist","; (fm^{#minus1}); #it{d}_{0} (fm)",50,-2,5);
    hist_emptyhist->GetXaxis()->SetRangeUser(-2, 5);
	hist_emptyhist->GetYaxis()->SetRangeUser(0.,18);

	hist_emptyhist->Draw();
	//pad1->SetLogz();
	hFinalEP_Rebinned->SetMaximum(0.5*Ndof);
	hFinalEP_Rebinned->SetMinimum(0*Ndof);
	hFinalEP_Rebinned->Draw("same,colz");
	hist_emptyhist->Draw("same,axis");

	//leg2->Draw("same");

    TLatex ScatLenText;
    ScatLenText.SetTextSize(hist_emptyhist->GetXaxis()->GetTitleSize());
    ScatLenText.SetNDC(kTRUE);
    ScatLenText.SetTextFont(12);
    ScatLenText.DrawLatex(0.810,0.042,"f_{#it{0}}^{#it{#minus1}}");
    ScatLenText.Draw();

    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.8);
    BeamText.SetNDC(kTRUE);

    BeamText.SetTextSize(gStyle->GetTextSize()*0.68);
    BeamText.DrawLatex(0.70, 0.36, "ALICE");
    BeamText.DrawLatex(0.70, 0.315, "pp #sqrt{#it{s}} = 7 TeV");
    //BeamText.DrawLatex(0.73, 0.315, "ALICE");
    BeamText.DrawLatex(0.70, 0.27, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.70, 0.225, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    BeamText.DrawLatex(0.70, 0.18, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs");
    BeamText.Draw();

	//pad1->RedrawAxis();
	//pad2->RedrawAxis();
	cLamLamExclusion->Write();


    Int_t MyPalette[100];
    Double_t Red[3]   = { 1,0,0 };
    Double_t Green[3]  = { 0,1,0 };
    Double_t Blue[3]    = { 0,0,1 };

    Double_t Length[2] = { 0, 1 };
    Int_t FI = TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, 100);
    for (int i=0;i<100;i++) MyPalette[i] = FI+i;


	//SetStyle3(false,0,NULL);
	gStyle->SetPalette(56);
	cLamLamExclusion->SaveAs(FolderName+"cLamLamChi2Exclusion.pdf");
	//SetStyle3(true,0,NULL);
	gStyle->SetPalette(56);
	cLamLamExclusion->SaveAs(FolderName+"cLamLamChi2Exclusion.png");

    delete [] nSigma_Contours;
    delete hist_emptyhist;
    delete grEntry1Fake;
    delete grEntry2Fake;
    delete grEntry3Fake;
    delete grEntry4Fake;
    delete hFinalEP_Rebinned;
    delete cLamLamExclusion;
    delete FileFinalEP;
    delete OutputFile;
}

void PlotDimiExclusionBE_ver1(const TString OutputFolder,
                              const TString FileNameStat, const TString FileNameSyst, const TString f0Inv_or_d0){

    if(f0Inv_or_d0!="f0Inv"&&f0Inv_or_d0!="d0"){
        printf("f0Inv or d0\nABORT!\n");
        return;
    }

    TFile* FileStat = new TFile(FileNameStat,"read");
    TH2F* hStat = (TH2F*)FileStat->Get(f0Inv_or_d0=="f0Inv"?"hEbin_f0Inv":"hEbin_d0");
    TGraph* grStatBest = (TGraph*)FileStat->Get(f0Inv_or_d0=="f0Inv"?"grBestFit_f0Inv":"grBestFit_d0");

    TFile* FileSyst = new TFile(FileNameSyst,"read");
    TH2F* hSyst = (TH2F*)FileSyst->Get(f0Inv_or_d0=="f0Inv"?"hEbin_f0Inv":"hEbin_d0");
    TGraph* grSystBest = (TGraph*)FileSyst->Get(f0Inv_or_d0=="f0Inv"?"grBestFit_f0Inv":"grBestFit_d0");

    TH2F* hSS = (TH2F*)hStat->Clone("hSS");
    for(unsigned uBinX=1; uBinX<=hSS->GetNbinsX(); uBinX++){
        for(unsigned uBinY=1; uBinY<=hSS->GetNbinsY(); uBinY++){
            if(hStat->GetBinContent(uBinX,uBinY)){
                hSS->SetBinContent(uBinX,uBinY,1);
            }
            else if(hSyst->GetBinContent(uBinX,uBinY)){
                hSS->SetBinContent(uBinX,uBinY,2);
            }
            else{
                hSS->SetBinContent(uBinX,uBinY,0);
            }
            hSS->SetBinError(uBinX,uBinY,0);
        }
    }
    const unsigned NumContours = 3;
    double* hSS_Contours = new double [NumContours];
    hSS_Contours[0] = 0;
    hSS_Contours[1] = 1.5;
    hSS_Contours[2] = 100;
    hSS->SetContour(NumContours,hSS_Contours);
    hStat->SetMarkerColor(fPalette2[1]);
    hSyst->SetMarkerColor(fPalette2[2]);
    hStat->SetFillColor(fPalette2[1]);
    hSyst->SetFillColor(fPalette2[2]);
    hStat->SetLineColor(fPalette2[1]);
    hSyst->SetLineColor(fPalette2[2]);

	const unsigned NumPlotObj = 3;
	TAttLine** ListOfObjects = new TAttLine* [NumPlotObj];
	ListOfObjects[0] = hStat;
	ListOfObjects[1] = hSyst;
	ListOfObjects[2] = hSS;
    SetStyle4(false,NumPlotObj,ListOfObjects);

    grStatBest->SetMarkerStyle(43);
    grStatBest->SetMarkerSize(4);
    grStatBest->SetMarkerColor(kRed-7);
    grStatBest->SetLineWidth(4);
    grStatBest->SetFillColor(kWhite);

    const float splitFraction = 0.8;
    TCanvas *cLamLamExclusionBE = new TCanvas("cLamLamExclusionBE", "cLamLamExclusionBE", 960, 540);
    TPad *pad1 = new TPad("lambda", "", 0.0,0.,splitFraction,1);
    TPad *pad2 = new TPad("lambda", "", splitFraction,0.,1,1);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    pad1->SetTopMargin(0.025);
    pad1->SetRightMargin(0.01);

    pad2->cd();
    TLegend* leg2 = new TLegend(0.01, 0.61, 0.99, 0.97);//lbrt
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(gStyle->GetTextSize()*2.25);
    leg2->AddEntry(hStat, TString::Format("#it{B}_{#Lambda#Lambda} (stat)"), "PF");
    leg2->AddEntry(hSyst, TString::Format("#it{B}_{#Lambda#Lambda} (syst)"), "PF");
    leg2->AddEntry(grStatBest, "Best fit", "P");
    leg2->Draw("same");

    pad1->cd();

    TH1F* hist_emptyhist;
    if(f0Inv_or_d0=="f0Inv"){
        hist_emptyhist = new TH1F("hist_emptyhist","; #it{B}_{#Lambda#Lambda} (MeV); #it{f}_{0}^{#minus1} (fm^{#minus1})",50,0,6.0);
        hist_emptyhist->GetYaxis()->SetRangeUser(0,6.0);
        hist_emptyhist->GetYaxis()->SetRangeUser(-0.6,0.);
    }
    else{
        hist_emptyhist = new TH1F("hist_emptyhist","; #it{B}_{#Lambda#Lambda} (MeV); #it{d}_{0} (fm)",50,0,6.0);
        hist_emptyhist->GetYaxis()->SetRangeUser(0,6.0);
        hist_emptyhist->GetYaxis()->SetRangeUser(0.,4.);
    }

	hist_emptyhist->Draw();
    hSS->Draw("same,col,hist");
    hist_emptyhist->Draw("same,axis");
    grStatBest->Draw("P,same");

    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.8);
    BeamText.SetNDC(kTRUE);

    double shift=0.5;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.68);
    BeamText.DrawLatex(0.75, 0.36+shift, "ALICE");
    BeamText.DrawLatex(0.75, 0.315+shift, "pp #sqrt{#it{s}} = 7 TeV");
    //BeamText.DrawLatex(0.73, 0.315, "ALICE");
    BeamText.DrawLatex(0.75, 0.27+shift, "pp #sqrt{#it{s}} = 13 TeV");
    BeamText.DrawLatex(0.75, 0.225+shift, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    BeamText.DrawLatex(0.75, 0.18+shift, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs");
    BeamText.Draw();

    //BeamText.Draw("same");
    if(f0Inv_or_d0=="f0Inv") cLamLamExclusionBE->SaveAs(OutputFolder+"cLamLam_BE_f0Inv.pdf");
    else cLamLamExclusionBE->SaveAs(OutputFolder+"cLamLam_BE_d0.pdf");

    delete hist_emptyhist;
    delete [] ListOfObjects;
    delete [] hSS_Contours;
    delete hSS;
    delete FileStat;
    delete FileSyst;
    delete leg2;
    delete pad1;
    delete pad2;
    delete cLamLamExclusionBE;

}

void plotCF(const TString system)
{

    char *expfile_pp = new char [256];
    char *expfile_pL = new char [256];
    char *expfile_LL = new char [256];
    char *expfile_pXim = new char [256];
    char *CATSfile = new char [256];
    char *simfile = new char [256];

    if(system=="pp"){
        strcpy(expfile_pp, "/home/dmihaylov/CernBox/Femto_pp13/data/AnalysisResults_AndiBernieSystME_CkProtonProton_0.root");
        strcpy(expfile_pL, "/home/dmihaylov/CernBox/Femto_pp13/data/AnalysisResults_AndiBernieSystME_CkProtonLambda_0.root");
        strcpy(expfile_LL, "/home/dmihaylov/CernBox/Femto_pp13/data/AnalysisResults_AndiBernieSystME_CkLambdaLambda_0.root");
        strcpy(expfile_pXim, "/home/dmihaylov/CernBox/Femto_pp13/data/AnalysisResults_AndiBernieSystME_CkProtonXim_0.root");
        strcpy(CATSfile, "/home/dmihaylov/Dropbox/FastForwardFemto/pp_13TeV/FitResults_110418/Gauss_LedniLo/SYSTEMATICS_CutVarAdd_Global_Radius_Normal.root");
    }
    else{
        strcpy(expfile_pp, "/home/dmihaylov/CernBox/pPb/DataFullest/AnalysisResults_AndiBernieSystME_CkProtonProton_0.root");
        strcpy(expfile_pL, "/home/dmihaylov/CernBox/pPb/DataFullest/AnalysisResults_AndiBernieSystME_CkProtonLambda_0.root");
        strcpy(expfile_LL, "/home/dmihaylov/CernBox/pPb/DataFullest/AnalysisResults_AndiBernieSystME_CkLambdaLambda_0.root");
        strcpy(expfile_pXim, "/home/dmihaylov/CernBox/pPb/DataFullest/AnalysisResults_AndiBernieSystME_CkProtonXim_0.root");
        strcpy(CATSfile, "/home/dmihaylov/Dropbox/FastForwardFemto/pp_13TeV/FitResults_110418/Gauss_LedniLo/SYSTEMATICS_CutVarAdd_Global_Radius_Normal.root");
    }
    strcpy(simfile, "");

  TDatime DateTime;
  TString Date = TString::Format("%i-%i-%i",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay());

  gStyle->SetCanvasPreferGL(1);
  const float right = 0.025;
  const float top = 0.025;

  // Mixed event normalisation
  const float normleft = 0.2;
  const float normright = 0.4;
  const float spinningDepth = 10.f;

  // for Data
  const float rebinData = 1;

  double energy;
  if(system=="pp"){
    energy=13;
  }
  else{
    energy=5.02;
  }

  bool EPOS = false;

  const char *prefix = "MB";
  const char *addon = "";

  TString data = "Data";
  TString sim = "DPMJET";
//  TString sim = "Pythia 8";

  //float r = 1.437;
  //float rErr = 0.011;
  //float rSystErrUp = 0.013;
  //float rSystErrDown = 0.006;
  float r = 1.185;
  float rErr = 0.009;
  float rSystErrUp = 0.016;
  float rSystErrDown = 0.009;

  float ppBL0 = 0.924;
  float ppBL1 = 0.386;
  float pLBL0 = 0.898;
  float pLBL1 = 0.355;
  float LLBL0 = 0.894;
  float LLBL1 = 0.320;
  float pXiBL0 = 0.919;
  float pXiBL1 = 0.453;

  //  EPOS
  if(EPOS) {
    r = 1.476;
    rErr = 0.015;
    rSystErrUp = 0.007;
    rSystErrDown = 0.014;

    ppBL0 = 0.931;
    ppBL1 = 0.399;
    pLBL0 = 0.938;
    pLBL1 = 0.328;
    LLBL0 = 0.893;
    LLBL1 = 0.321;
    pXiBL0 = 0.921;
    pXiBL1 = 0.443;
  }
  SetStyle();

    TFile* file_pp = new TFile(expfile_pp, "read");
    TH1F* hCk_pp = (TH1F*)file_pp->Get("hCkTotNormWeight");

    TFile* file_pL = new TFile(expfile_pL, "read");
    TH1F* hCk_pL = (TH1F*)file_pL->Get("hCkTotNormWeight");

    TFile* file_LL = new TFile(expfile_LL, "read");
    TH1F* hCk_LL = (TH1F*)file_LL->Get("hCkTotNormWeight");

    TFile* file_pXim = new TFile(expfile_pXim, "read");
    TH1F* hCk_pXim = (TH1F*)file_pXim->Get("hCkTotNormWeight");
    //hCk_ReweightedMeV_1

  // EXP DATA
  TFile* _file0=TFile::Open("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample6/AnalysisResults.root");
  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon),Results);
  TList* tmpFolder=(TList*)Results->FindObject("Particle0_Particle0");
  TH1F* histRE_relK_pp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");
//  histRE_relK_pp->Rebin(rebinData);
  TH1F* histME_relK_pp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");
//  histME_relK_pp->Rebin(rebinData);
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle1");
  TH1F* histRE_relK_ApAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle1");
//  histRE_relK_ApAp->Rebin(rebinData);
  TH1F* histME_relK_ApAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle1");
//  histME_relK_ApAp->Rebin(rebinData);
  tmpFolder=(TList*)Results->FindObject("Particle0_Particle2");
  TH1F* histRE_relK_Lp = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");
  histRE_relK_Lp->Rebin(rebinData);
  TH1F* histME_relK_Lp = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");
  histME_relK_Lp->Rebin(rebinData);
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle3");
  TH1F* histRE_relK_ALAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle3");
  histRE_relK_ALAp->Rebin(rebinData);
  TH1F* histME_relK_ALAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle3");
  histME_relK_ALAp->Rebin(rebinData);
  tmpFolder=(TList*)Results->FindObject("Particle2_Particle2");
  TH1F* histRE_relK_LL = (TH1F*)tmpFolder->FindObject("SEDist_Particle2_Particle2");
  histRE_relK_LL->Rebin(rebinData);
  TH1F* histME_relK_LL = (TH1F*)tmpFolder->FindObject("MEDist_Particle2_Particle2");
  histME_relK_LL->Rebin(rebinData);
  tmpFolder=(TList*)Results->FindObject("Particle3_Particle3");
  TH1F* histRE_relK_ALAL = (TH1F*)tmpFolder->FindObject("SEDist_Particle3_Particle3");
  histRE_relK_ALAL->Rebin(rebinData);
  TH1F* histME_relK_ALAL = (TH1F*)tmpFolder->FindObject("MEDist_Particle3_Particle3");
  histME_relK_ALAL->Rebin(rebinData);
  tmpFolder=(TList*)Results->FindObject("Particle0_Particle4");
  TH1F* histRE_relK_Xip = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle4");
  histRE_relK_Xip->Rebin(rebinData);
  TH1F* histME_relK_Xip = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle4");
  histME_relK_Xip->Rebin(rebinData);
  tmpFolder=(TList*)Results->FindObject("Particle1_Particle5");
  TH1F* histRE_relK_AXiAp = (TH1F*)tmpFolder->FindObject("SEDist_Particle1_Particle5");
  histRE_relK_AXiAp->Rebin(rebinData);
  TH1F* histME_relK_AXiAp = (TH1F*)tmpFolder->FindObject("MEDist_Particle1_Particle5");
  histME_relK_AXiAp->Rebin(rebinData);
  TH1F *hist_CF_Lp_ALAp_exp[3];
  TH1F *hist_CF_LL_ALAL_exp[3];
  TH1F *hist_CF_pp_ApAp_exp[3];
  TH1F *hist_CF_pXi_ApAXi_exp[3];

  hist_CF_Lp_ALAp_exp[0] = Calculate_CF(histRE_relK_Lp,histME_relK_Lp,"hist_CF_Lp_exp",normleft,normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[1] = Calculate_CF(histRE_relK_ALAp,histME_relK_ALAp,"hist_CF_ALAp_exp",normleft,normright, addon, spinningDepth);
  hist_CF_Lp_ALAp_exp[2] = add_CF(hist_CF_Lp_ALAp_exp[0],hist_CF_Lp_ALAp_exp[1],"hist_CF_Lp_ALAp_exp_sum");
  hist_CF_LL_ALAL_exp[0] = Calculate_CF(histRE_relK_LL,histME_relK_LL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[1] = Calculate_CF(histRE_relK_ALAL,histME_relK_ALAL,"hist_CF_LL_exp",normleft,normright, addon, spinningDepth);
  hist_CF_LL_ALAL_exp[2] = add_CF(hist_CF_LL_ALAL_exp[0],hist_CF_LL_ALAL_exp[1],"hist_CF_LL_ALAL_exp_sum");
  hist_CF_pp_ApAp_exp[0] = Calculate_CF(histRE_relK_pp,histME_relK_pp,"hist_CF_pp",normleft,normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[1] = Calculate_CF(histRE_relK_ApAp,histME_relK_ApAp,"hist_CF_ApAp",normleft,normright, addon, spinningDepth);
  hist_CF_pp_ApAp_exp[2] = add_CF(hist_CF_pp_ApAp_exp[0],hist_CF_pp_ApAp_exp[1],"hist_CF_pp_ApAp_exp_sum");
  hist_CF_pXi_ApAXi_exp[0] = Calculate_CF(histRE_relK_Xip,histME_relK_Xip,"hist_CF_pXi",normleft,normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[1] = Calculate_CF(histRE_relK_AXiAp,histME_relK_AXiAp,"hist_CF_ApAXi",normleft,normright, addon, spinningDepth);
  hist_CF_pXi_ApAXi_exp[2] = add_CF(hist_CF_pXi_ApAXi_exp[0],hist_CF_pXi_ApAXi_exp[1],"hist_CF_pXi_ApAXi_exp_sum");

    for(unsigned uBin=0; uBin<60; uBin++){
        hist_CF_pp_ApAp_exp[2]->SetBinContent(uBin+1,hCk_pp->GetBinContent(uBin+1));
        hist_CF_pp_ApAp_exp[2]->SetBinError(uBin+1,hCk_pp->GetBinError(uBin+1));
    }
    hCk_pp = hist_CF_pp_ApAp_exp[2];

    for(unsigned uBin=0; uBin<20; uBin++){
        hist_CF_Lp_ALAp_exp[2]->SetBinContent(uBin+1,hCk_pL->GetBinContent(uBin+1));
        hist_CF_Lp_ALAp_exp[2]->SetBinError(uBin+1,hCk_pL->GetBinError(uBin+1));
    }
    hCk_pL = hist_CF_Lp_ALAp_exp[2];

    for(unsigned uBin=0; uBin<60; uBin++){
        hist_CF_pp_ApAp_exp[2]->SetBinContent(uBin+1,hCk_pp->GetBinContent(uBin+1));
        hist_CF_pp_ApAp_exp[2]->SetBinError(uBin+1,hCk_pp->GetBinError(uBin+1));
    }
    hCk_pp = hist_CF_pp_ApAp_exp[2];

    for(unsigned uBin=0; uBin<20; uBin++){
        hist_CF_LL_ALAL_exp[2]->SetBinContent(uBin+1,hCk_LL->GetBinContent(uBin+1));
        hist_CF_LL_ALAL_exp[2]->SetBinError(uBin+1,hCk_LL->GetBinError(uBin+1));
    }
    hCk_LL = hist_CF_LL_ALAL_exp[2];

    for(unsigned uBin=0; uBin<20; uBin++){
        hist_CF_pXi_ApAXi_exp[2]->SetBinContent(uBin+1,hCk_pXim->GetBinContent(uBin+1));
        hist_CF_pXi_ApAXi_exp[2]->SetBinError(uBin+1,hCk_pXim->GetBinError(uBin+1));
    }
    hCk_pXim = hist_CF_pXi_ApAXi_exp[2];

  SetStyleHisto(hCk_pp, 0,0);
  SetStyleHisto(hCk_pL, 0,0);
  SetStyleHisto(hCk_LL, 0,0);
  SetStyleHisto(hCk_pXim, 0,0);

  // SYSTEMATIC UNCERTAINTIES
  TString input_sys_pp;
  TString input_sys_pL;
  TString input_sys_LL;
  TString input_sys_pXi;
  if(system=="pp"){
      input_sys_pp = "/home/dmihaylov/CernBox/Femto_pp13/C2totalsysPP.root";
      input_sys_pL = "/home/dmihaylov/CernBox/Femto_pp13/C2totalsysPL.root";
      input_sys_LL = "/home/dmihaylov/CernBox/Femto_pp13/C2totalsysLL.root";
      input_sys_pXi = "/home/dmihaylov/CernBox/Femto_pp13/C2totalsysPXi.root";
  }
  else{
      input_sys_pp = "";
      input_sys_pL = "";
      input_sys_LL = "";
      input_sys_pXi = "";
  }
  TFile* file_sys_pp = new TFile(input_sys_pp.Data());
  TFile* file_sys_pL = new TFile(input_sys_pL.Data());
  TFile* file_sys_pXi = new TFile(input_sys_pXi.Data());
  TFile* file_sys_LL = new TFile(input_sys_LL.Data());
  TH1F* hist_sys_pp = (TH1F*)file_sys_pp->Get("C2totalsysPP");
  hist_sys_pp->SetLineWidth(2.0);
  TH1F* hist_sys_pL = (TH1F*)file_sys_pL->Get("C2totalsysPL");
  hist_sys_pL->SetLineWidth(2.0);
  TH1F* hist_sys_LL = (TH1F*)file_sys_LL->Get("C2totalsysLL");
  hist_sys_LL->SetLineWidth(2.0);
  TH1F* hist_sys_pXi = (TH1F*)file_sys_pXi->Get("C2totalsysPXi");
  hist_sys_pXi->SetLineWidth(2.0);

  // UNCERTAINTIES ON THE FIT
  TFile *systFit = TFile::Open(CATSfile);
  TGraph *grppDefault = nullptr;
  TGraph *grppLow = nullptr;
  TGraph *grppUp = nullptr;
  TGraph *grpLDefaultNLO = nullptr;
  TGraph *grpLLowNLO = nullptr;
  TGraph *grpLUpNLO = nullptr;
  TGraph *grpLDefaultLO = nullptr;
  TGraph *grpLLowLO = nullptr;
  TGraph *grpLUpLO = nullptr;
  TGraph *grLLDefault = nullptr;
  TGraph *grLLLow = nullptr;
  TGraph *grLLUp = nullptr;
  TGraph *grpXiDefault = nullptr;
  TGraph *grpXiLower = nullptr;
  TGraph *grpXiUpper = nullptr;
  TGraph *grpXiDefaultCoulomb = nullptr;
  TGraph *grpXiLowerCoulomb = nullptr;
  TGraph *grpXiUpperCoulomb = nullptr;
  TGraphErrors *grFemtopp = nullptr;
  TGraphErrors *grFemtopLNLO = nullptr;
  TGraphErrors *grFemtopLLO = nullptr;
  TGraphErrors *grFemtoLL = nullptr;
  TGraphErrors *grFemtopXi = nullptr;
  TGraphErrors *grFemtopXiCoulomb = nullptr;
  if(systFit) {
    grppDefault = (TGraph*)systFit->Get("ppGraphDefault");
    grppLow = (TGraph*)systFit->Get("ppGraphLowerLim");
    grppUp = (TGraph*)systFit->Get("ppGraphUpperLim");

    const char* nlostring = (EPOS) ? "" : "_NLO";
    const char* lostring = (EPOS) ? "_LO" : "";
    const char* prefix = (EPOS) ? "" : "Copy_";
    const char* prefixLO = (!EPOS) ? "" : "Copy_";
    grpLDefaultNLO = (TGraph*)systFit->Get(Form("pLamGraphDefault%s", nlostring));
    grpLLowNLO = (TGraph*)systFit->Get(Form("%spLamGraphLowerLim%s", prefix, nlostring));
    grpLUpNLO = (TGraph*)systFit->Get(Form("pLamGraphUpperLim%s", nlostring));
    grpLDefaultLO = (TGraph*)systFit->Get(Form("pLamGraphDefault%s", lostring));
    grpLLowLO = (TGraph*)systFit->Get(Form("%spLamGraphLowerLim%s", prefixLO, lostring));
    grpLUpLO = (TGraph*)systFit->Get(Form("pLamGraphUpperLim%s", lostring));

    grLLDefault = (TGraph*)systFit->Get("LamLamGraphDefault");
    grLLLow = (TGraph*)systFit->Get("LamLamGraphLowerLim");
    grLLUp = (TGraph*)systFit->Get("LamLamGraphUpperLim");
    grpXiDefault = (TGraph*)systFit->Get("pXimGraphDefault");
    grpXiLower = (TGraph*)systFit->Get("pXimGraphLowerLim");
    grpXiUpper = (TGraph*)systFit->Get("pXimGraphUpperLim");

    grpXiDefaultCoulomb = (TGraph*)systFit->Get("pXimGraphDefault_COULOMB");
    grpXiLowerCoulomb = (TGraph*)systFit->Get("pXimGraphLowerLim_COULOMB");
    grpXiUpperCoulomb = (TGraph*)systFit->Get("pXimGraphUpperLim_COULOMB");

    grFemtopp = FemtoModelFitBands(grppDefault, grppLow, grppUp);
    grFemtopp->SetFillColor(fColors[2]);
    grFemtopp->SetLineColor(fColors[2]);
    grFemtopp->SetLineWidth(3);
    grFemtopLNLO = FemtoModelFitBands(grpLDefaultNLO, grpLLowNLO, grpLUpNLO);
    grFemtopLNLO->SetFillColor(fColors[1]);
    grFemtopLNLO->SetLineColor(fColors[1]);
    grFemtopLNLO->SetLineWidth(3);
    grFemtopLLO = FemtoModelFitBands(grpLDefaultLO, grpLLowLO, grpLUpLO);
    grFemtopLLO->SetFillColor(fColors[3]);
    grFemtopLLO->SetLineColor(fColors[3]);
    grFemtopLLO->SetLineWidth(3);
    grFemtoLL = FemtoModelFitBands(grLLDefault, grLLLow, grLLUp);
    grFemtoLL->SetFillColor(fColors[5]);
    grFemtoLL->SetLineColor(fColors[5]);
    grFemtoLL->SetLineWidth(3);
    grFemtopXi = FemtoModelFitBands(grpXiDefault, grpXiLower, grpXiUpper);
    grFemtopXi->SetFillColor(fColors[6]);
    grFemtopXi->SetLineColor(fColors[6]);
    grFemtopXi->SetLineWidth(3);
    grFemtopXiCoulomb = FemtoModelFitBands(grpXiDefaultCoulomb, grpXiLowerCoulomb, grpXiUpperCoulomb);
    grFemtopXiCoulomb->SetFillColor(fColors[7]);
    grFemtopXiCoulomb->SetLineColor(fColors[7]);
    grFemtopXiCoulomb->SetLineWidth(3);
  }

  TGraph *grFakeSyst = new TGraph();
  SetStyleGraph(grFakeSyst,0,0);
  grFakeSyst->SetFillColor(fFillColors[0]);
  grFakeSyst->SetLineColor(fFillColors[0]);
  TGraph *grFakePP = new TGraph();
  grFakePP->SetLineColor(fColors[2]);
  grFakePP->SetLineWidth(4);
  TGraph *grFakePLnlo = new TGraph();
  grFakePLnlo->SetLineColor(fColors[1]);
  grFakePLnlo->SetLineWidth(4);
  TGraph *grFakePLlo = new TGraph();
  grFakePLlo->SetLineColor(fColors[3]);
  grFakePLlo->SetLineWidth(4);
  TGraph *grFakeLL = new TGraph();
  grFakeLL->SetLineColor(fColors[5]);
  grFakeLL->SetLineWidth(4);
  TGraph *grFakeLLStar = new TGraph();
  grFakeLLStar->SetLineColor(fColors[6]);
  grFakeLLStar->SetLineWidth(4);
  TGraph *grFakepXi= new TGraph();
  grFakepXi->SetLineColor(fColors[6]);
  grFakepXi->SetLineWidth(4);
  TGraph *grFakeXiCoulomb= new TGraph();
  grFakeXiCoulomb->SetLineColor(fColors[7]);
  grFakeXiCoulomb->SetLineWidth(4);

  // BASELINE
  TF1 *baselinePP = new TF1("baselinePP", "pol1", 0, 1);
  TF1 *baselinePL = new TF1("baselinePL", "pol1", 0, 1);
  TF1 *baselineLL = new TF1("baselineLL", "pol1", 0, 1);
  TF1 *baselinePXI = new TF1("baselinePXI", "pol1", 0, 1);
  baselinePP->SetParameter(0, ppBL0);
  baselinePP->SetParameter(1, ppBL1);
  baselinePL->SetParameter(0, pLBL0);
  baselinePL->SetParameter(1, pLBL1);
  baselineLL->SetParameter(0, LLBL0);
  baselineLL->SetParameter(1, LLBL1);
  baselinePXI->SetParameter(0, pXiBL0);
  baselinePXI->SetParameter(1, pXiBL1);
  baselinePP->SetLineStyle(2);
  baselinePL->SetLineStyle(2);
  baselineLL->SetLineStyle(2);
  baselinePXI->SetLineStyle(2);
  baselinePP->SetLineColor(fFillColors[6]);
  baselinePL->SetLineColor(fFillColors[6]);
  baselineLL->SetLineColor(fFillColors[6]);
  baselinePXI->SetLineColor(fFillColors[6]);

  TF1* qsLL = new TF1("qsLL","([0]+[1]*x)*(1.-[2]*0.5*exp(-4.*x*x*[3]*[3]*5.0677*5.0677))", 0, 1);
  qsLL->SetParameter(0, LLBL0);
  qsLL->SetParameter(1, LLBL1);
  qsLL->SetParameter(2, 0.36);//THE LAMBDA PARAMETER!!!
  qsLL->SetParameter(3, r);
  qsLL->SetLineColor(kGray+1);
  qsLL->SetLineStyle(5);
  qsLL->SetLineWidth(3);

  // PLOTTING

  TCanvas *Can_CF_fitting = new TCanvas("Can_CF_fitting","Can_CF_fitting",0,0,1100,1000);
  Can_CF_fitting->Divide(2,2);
  Can_CF_fitting->cd(1);
  Can_CF_fitting->cd(1)->SetRightMargin(right);
  Can_CF_fitting->cd(1)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_pp_ApAp = DrawSystematicError(hCk_pp, hist_sys_pp, 0.002);
  Tgraph_syserror_pp_ApAp->SetLineColor(kWhite);
  Tgraph_syserror_pp_ApAp->Draw("Ap");
  baselinePP->Draw("same");
  Tgraph_syserror_pp_ApAp->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pp_ApAp->GetXaxis()->SetRangeUser(0, 0.125);
  Tgraph_syserror_pp_ApAp->GetYaxis()->SetRangeUser(0.5, 3.5);
  if(grFemtopp) grFemtopp->Draw("L3 same");
  Tgraph_syserror_pp_ApAp->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pp_ApAp->Draw("2 same");
  hCk_pp->Draw("pe same");
  TLegend *legpp = new TLegend(0.48,0.595,0.85,0.81);
  legpp->SetBorderSize(0);
  legpp->SetTextFont(42);
  legpp->SetTextSize(gStyle->GetTextSize()*0.75);
  legpp->AddEntry(hCk_pp, "pp #oplus #bar{p}#bar{p} pairs", "pe");
  legpp->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legpp->AddEntry(grFakePP,"Femtoscopic fit","l");
  legpp->Draw("same");
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);
  if(system=="pp") BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  TLatex text;
  text.SetTextSize(gStyle->GetTextSize()*0.75);
  text.SetNDC();
  text.SetTextColor(1);
  if(!EPOS) text.DrawLatex(0.5, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.5, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));

  Can_CF_fitting->cd(2);
  Can_CF_fitting->cd(2)->SetRightMargin(right);
  Can_CF_fitting->cd(2)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_pL_ApAL = DrawSystematicError(hCk_pL, hist_sys_pL, 0.005);
  Tgraph_syserror_pL_ApAL->SetLineColor(kWhite);
  Tgraph_syserror_pL_ApAL->Draw("ap");
  baselinePL->Draw("same");
  Tgraph_syserror_pL_ApAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pL_ApAL->GetYaxis()->SetRangeUser(0.8, 2.);
  if(grFemtopLNLO) grFemtopLNLO->Draw("l3 same");
  if(!EPOS && grFemtopLLO) grFemtopLLO->Draw("l3 same");
  Tgraph_syserror_pL_ApAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pL_ApAL->Draw("2 same");
  hCk_pL->Draw("pe same");
  TLegend *legLp;
  if(EPOS) legLp = new TLegend(0.385,0.595,0.75,0.81);
  else legLp     = new TLegend(0.385,0.545,0.75,0.81);
  legLp->SetBorderSize(0);
  legLp->SetTextFont(42);
  legLp->SetTextSize(gStyle->GetTextSize()*0.75);
  legLp->AddEntry(hCk_pL, "p#Lambda #oplus #bar{p}#bar{#Lambda} pairs", "pe");
  legLp->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legLp->AddEntry(grFakePLnlo,"Femtoscopic fit (#chiEFT NLO)","l");
  if(!EPOS) legLp->AddEntry(grFakePLlo,"Femtoscopic fit (#chiEFT LO)","l");
  legLp->Draw("same");
  if(system=="pp") BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  if(!EPOS) text.DrawLatex(0.4, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.4, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize()*0.6);
  ref.SetNDC(kTRUE);
  if(EPOS) ref.DrawLatex(0.48, 0.565, "Nucl. Phys. A915 (2013) 24.");
  else ref.DrawLatex(0.48, 0.515, "Nucl. Phys. A915 (2013) 24.");

  Can_CF_fitting->cd(3);
  Can_CF_fitting->cd(3)->SetRightMargin(right);
  Can_CF_fitting->cd(3)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError(hCk_LL, hist_sys_LL, 0.005);
  Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);
  Tgraph_syserror_LL_ALAL->Draw("ap");
  baselineLL->Draw("same");
  qsLL->Draw("same");
  Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.35, 2.);
  if(grFemtoLL) grFemtoLL->Draw("l3 same");
  Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_LL_ALAL->Draw("2 same");
  hCk_LL->Draw("pe same");
  TLegend *legLL = new TLegend(0.48,0.595,0.85,0.81);
  legLL->SetBorderSize(0);
  legLL->SetTextFont(42);
  legLL->SetTextSize(gStyle->GetTextSize()*0.75);
  legLL->AddEntry(hCk_pp, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs", "pe");
  legLL->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legLL->AddEntry(grFakeLL,"Femtoscopic fit","l");
  legLL->AddEntry(qsLL,"Quantum statistics","l");
  legLL->Draw("same");
  if(system=="pp") BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.5, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  if(!EPOS) text.DrawLatex(0.5, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.5, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));

  Can_CF_fitting->cd(4);
  Can_CF_fitting->cd(4)->SetRightMargin(right);
  Can_CF_fitting->cd(4)->SetTopMargin(top);
  TGraphErrors *Tgraph_syserror_pXi_ApAXi = DrawSystematicError(hCk_pXim, hist_sys_pXi, 0.005);
  Tgraph_syserror_pXi_ApAXi->SetLineColor(kWhite);
  Tgraph_syserror_pXi_ApAXi->Draw("ap");
  baselinePXI->Draw("same");
  Tgraph_syserror_pXi_ApAXi->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pXi_ApAXi->GetYaxis()->SetRangeUser(0.75, 3.);
  Tgraph_syserror_pXi_ApAXi->SetFillColorAlpha(kBlack, 0.4);
  if(grFemtopXi) grFemtopXi->Draw("l3 same");
  if(grFemtopXiCoulomb) grFemtopXiCoulomb->Draw("l3 same");
  Tgraph_syserror_pXi_ApAXi->Draw("2 same");
  hCk_pXim->Draw("pe same");
  TLegend *legpXi = new TLegend(0.385,0.545,0.75,0.81);
  legpXi->SetBorderSize(0);
  legpXi->SetTextFont(42);
  legpXi->SetTextSize(gStyle->GetTextSize()*0.75);
  legpXi->AddEntry(hCk_pp, "p#Xi^{-} #oplus #bar{p}#Xi^{+} pairs", "pe");
  legpXi->AddEntry(grFakeSyst, "Syst. uncertainties", "f");
  legpXi->AddEntry(grFakeXiCoulomb, "Femtoscopic fit (Coulomb)","l");
  legpXi->AddEntry(grFakepXi,"Femtoscopic fit (HAL QCD)","l");
  legpXi->Draw("same");
  ref.DrawLatex(0.48, 0.515, "Nucl. Phys. A967 (2017) 856.");
  ref.DrawLatex(0.48, 0.475, "PoS LATTICE2016 (2017) 116.");
  if(system=="pp") BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.4, 0.875, Form("ALICE %s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  if(!EPOS) text.DrawLatex(0.4, 0.825, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.4, 0.825, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_fitting->Print(TString::Format("./FemtoBoyzPlots/%sCF_EPOS_prelim.pdf",Date.Data()));
  else Can_CF_fitting->Print(TString::Format("./FemtoBoyzPlots/%s-CF_Gauss_prelim.pdf",Date.Data()));

  // PRELIMINARY PLOTS
  //gStyle->SetErrorX(0.);
  gStyle->SetErrorX(0.001);
  TCanvas *Can_CF_pp = new TCanvas("pp","pp", 0,0,650,550);
  Can_CF_pp->SetRightMargin(right);
  Can_CF_pp->SetTopMargin(top);
  Tgraph_syserror_pp_ApAp->SetLineColor(kWhite);
  Tgraph_syserror_pp_ApAp->Draw("Ap");
  baselinePP->Draw("same");
  Tgraph_syserror_pp_ApAp->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pp_ApAp->GetXaxis()->SetRangeUser(0, 0.125);
  Tgraph_syserror_pp_ApAp->GetYaxis()->SetRangeUser(0.5, 3.5);
  if(grFemtopp) grFemtopp->Draw("L3 same");
  Tgraph_syserror_pp_ApAp->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pp_ApAp->Draw("2 same");
  hCk_pp->Draw("pe same");
  TLegend *legpp2 = new TLegend(0.53,0.545,0.85,0.76);
  legpp2->SetBorderSize(0);
  legpp2->SetTextFont(42);
  legpp2->SetTextSize(gStyle->GetTextSize()*0.75);
  legpp2->AddEntry(grFakeSyst, "pp #oplus #bar{p}#bar{p} pairs", "fpe");
//  legpp2->AddEntry(hist_CF_pp_ApAp_exp[2], "with Syst. uncertainties", "");
  legpp2->AddEntry(baselinePP,"Baseline","l");
  legpp2->AddEntry(grFakePP,"Femtoscopic fit","l");
  legpp2->Draw("same");
  BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.55, 0.875, "ALICE");
  if(system=="pp") BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  text.SetTextSize(gStyle->GetTextSize()*0.75);
  text.SetNDC();
  text.SetTextColor(1);
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_pp->Print(TString::Format("./FemtoBoyzPlots/%s-CF_pp_EPOS_prelim.pdf",Date.Data()));
  else Can_CF_pp->Print(TString::Format("./FemtoBoyzPlots/%s-CF_pp_Gauss_prelim.pdf",Date.Data()));


  TCanvas *Can_CF_pL = new TCanvas("pL","pL", 0,0,650,550);
  Can_CF_pL->SetRightMargin(right);
  Can_CF_pL->SetTopMargin(top);
  Tgraph_syserror_pL_ApAL->SetLineColor(kWhite);
  Tgraph_syserror_pL_ApAL->Draw("ap");
  baselinePL->Draw("same");
  Tgraph_syserror_pL_ApAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pL_ApAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pL_ApAL->GetYaxis()->SetRangeUser(0.8, 2.);
  if(grFemtopLNLO) grFemtopLNLO->Draw("l3 same");
  if(!EPOS && grFemtopLLO) grFemtopLLO->Draw("l3 same");
  Tgraph_syserror_pL_ApAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_pL_ApAL->Draw("2 same");
  hCk_pL->Draw("pe same");
  TLegend *legLp2;
  if(EPOS) legLp2 = new TLegend(0.53,0.545,0.85,0.76);
  else legLp2 = new TLegend(0.53,0.495,0.85,0.76);
  legLp2->SetBorderSize(0);
  legLp2->SetTextFont(42);
  legLp2->SetTextSize(gStyle->GetTextSize()*0.75);
  legLp2->AddEntry(grFakeSyst, "p#Lambda #oplus #bar{p}#bar{#Lambda} pairs", "fpe");
//  legLp2->AddEntry(hist_CF_Lp_ALAp_exp[2], "with Syst. uncertainties", "");
  legLp2->AddEntry(baselinePL,"Baseline","l");
  legLp2->AddEntry(grFakePLnlo,"Femtoscopic fit (#chiEFT NLO)","l");
  if(!EPOS) legLp2->AddEntry(grFakePLlo,"Femtoscopic fit (#chiEFT LO)","l");
  legLp2->Draw("same");
  BeamText.DrawLatex(0.55, 0.875, "ALICE");
  if(system=="pp") BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) ref.DrawLatex(0.61, 0.515, "Nucl. Phys. A915 (2013) 24");
  else ref.DrawLatex(0.61, 0.465, "Nucl. Phys. A915 (2013) 24");
  if(EPOS) Can_CF_pL->Print(TString::Format("./FemtoBoyzPlots/%s-CF_pL_EPOS_prelim.pdf",Date.Data()));
  else Can_CF_pL->Print(TString::Format("./FemtoBoyzPlots/%s-CF_pL_Gauss_prelim.pdf",Date.Data()));

  TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
  Can_CF_LL->SetRightMargin(right);
  Can_CF_LL->SetTopMargin(top);
  Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);
  Tgraph_syserror_LL_ALAL->Draw("ap");
  baselineLL->Draw("same");
  qsLL->Draw("same");
  Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.35, 2.);
  if(grFemtoLL) grFemtoLL->Draw("l3 same");
  Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror_LL_ALAL->Draw("2 same");
  hCk_LL->Draw("pe same");
  TLegend *legLL2 = new TLegend(0.53,0.545,0.85,0.76);
  legLL2->SetBorderSize(0);
  legLL2->SetTextFont(42);
  legLL2->SetTextSize(gStyle->GetTextSize()*0.75);
  legLL2->AddEntry(grFakeSyst, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs", "fpe");
//  legLL2->AddEntry(hist_CF_LL_ALAL_exp[2], "with Syst. uncertainties", "");
  legLL2->AddEntry(baselineLL,"Baseline","l");
  legLL2->AddEntry(grFakeLL,"Femtoscopic fit","l");
  legLL2->AddEntry(qsLL,"Quantum statistics","l");
  legLL2->Draw("same");
  BeamText.DrawLatex(0.55, 0.875, "ALICE");
  if(system=="pp") BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_LL->Print(TString::Format("./FemtoBoyzPlots/%s-CF_LL_EPOS_prelim.pdf",Date.Data()));
  else Can_CF_LL->Print(TString::Format("./FemtoBoyzPlots/%s-CF_LL_Gauss_prelim.pdf",Date.Data()));

  TCanvas *Can_CF_pXi = new TCanvas("pXi","pXi", 0,0,650,550);
  Can_CF_pXi->SetRightMargin(right);
  Can_CF_pXi->SetTopMargin(top);
  Tgraph_syserror_pXi_ApAXi->SetLineColor(kWhite);
  Tgraph_syserror_pXi_ApAXi->Draw("ap");
  baselinePXI->Draw("same");
  Tgraph_syserror_pXi_ApAXi->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetRangeUser(0, 0.2);
  Tgraph_syserror_pXi_ApAXi->GetXaxis()->SetNdivisions(505);
  Tgraph_syserror_pXi_ApAXi->GetYaxis()->SetRangeUser(0.75, 5.);
  Tgraph_syserror_pXi_ApAXi->SetFillColorAlpha(kBlack, 0.4);
  if(grFemtopXi) grFemtopXi->Draw("l3 same");
  if(grFemtopXiCoulomb) grFemtopXiCoulomb->Draw("l3 same");
  Tgraph_syserror_pXi_ApAXi->Draw("2 same");
  hCk_pXim->Draw("pe same");
  TLegend *legpXi2 = new TLegend(0.53,0.495,0.85,0.76);
  legpXi2->SetBorderSize(0);
  legpXi2->SetTextFont(42);
  legpXi2->SetTextSize(gStyle->GetTextSize()*0.75);
  legpXi2->AddEntry(grFakeSyst, "p#Xi^{-} #oplus #bar{p}#Xi^{+} pairs", "fpe");
//  legpXi2->AddEntry(hist_CF_pXi_ApAXi_exp[2], "with Syst. uncertainties", "");
  legpXi2->AddEntry(baselinePXI,"Baseline","l");
  legpXi2->AddEntry(grFakeXiCoulomb, "Femtoscopic fit (Coulomb)","l");
  legpXi2->AddEntry(grFakepXi,"Femtoscopic fit (HAL QCD)","l");
  legpXi2->Draw("same");

  ref.DrawLatex(0.61, 0.465, "Nucl. Phys. A967 (2017) 856");
  ref.DrawLatex(0.61, 0.415, "PoS LATTICE2016 (2017) 116");
  BeamText.DrawLatex(0.55, 0.875, "ALICE");
  if(system=="pp") BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
  else BeamText.DrawLatex(0.55, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", system.Data(), energy));
  if(!EPOS) text.DrawLatex(0.55, 0.775, Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", r, rErr, rSystErrUp, rSystErrDown));
  else text.DrawLatex(0.55, 0.775, Form("#it{N}_{R} = %.3f #pm %.3f ^{+%.3f}_{-%.3f}", r, rErr, rSystErrUp, rSystErrDown));
  if(EPOS) Can_CF_pXi->Print(TString::Format("./FemtoBoyzPlots/%s-CF_pXi_EPOS_prelim.pdf",Date.Data()));
  else Can_CF_pXi->Print(TString::Format("./FemtoBoyzPlots/%s-CF_pXi_Gauss_prelim.pdf",Date.Data()));

  auto *outfile = new TFile("./FemtoBoyzPlots/CFout.root", "RECREATE");
  hCk_pp->Write("pp");
  hCk_pL->Write("pL");
  hCk_LL->Write("LL");
  hCk_pXim->Write("pXi");

}

//Data set 0: a=9.584e-01; b=1.247e-04
//Data set 1: a=1.001e+00; b=-1.067e-05
//Data set 2: a=1.003e+00; b=-2.675e-06
void PlotLamLamFit(const TString FitFolder, const unsigned WhichDataSet, const unsigned NumSystIter,
                   const TString DataFileName, const TString DataHistoName, const double parA, const double parB,
                   const TString DataSystFileName){
    TFile* DataFile = new TFile(DataFileName,"read");
    TH1F* DataHisto = (TH1F*)DataFile->Get(DataHistoName);
    TFile* DataSystFile = new TFile(DataSystFileName);
    TH1F* DataSystHisto = (TH1F*)DataSystFile->Get("C2totalsysLL");

    TGraph gDownLimit;
    TGraph gUpLimit;

    gStyle->SetCanvasPreferGL(1);
    SetStyle();
    const float right = 0.025;
    const float top = 0.025;

    for(unsigned uSyst=0; uSyst<NumSystIter; uSyst++){
//printf("hi\n");
        TFile* SystFile = new TFile(FitFolder+TString::Format("SystId%u_SystJobId%u_JobId0_OutputLamLamConfidence.root",uSyst,uSyst),"read");
//printf("SystFile=%p\n",SystFile);
//SystFile->ls();
        TGraph* gFitResult = (TGraph*)SystFile->Get(TString::Format("grFitRes_%u_%u",WhichDataSet,uSyst));
//printf("gFitResult=%p (%s)\n",gFitResult,TString::Format("grFitRes_%u_%u",WhichDataSet).Data());
        unsigned NumBins = gFitResult->GetN();
        double MiniumValue;
        double MaximumValue;
        double CurrentValue;
        double Momentum;
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            if(uSyst==0){
                gFitResult->GetPoint(uBin,Momentum,CurrentValue);
                gDownLimit.SetPoint(uBin,Momentum,CurrentValue);

                gUpLimit.SetPoint(uBin,Momentum,CurrentValue);

            }
            else{
                gFitResult->GetPoint(uBin,Momentum,CurrentValue);
                gDownLimit.GetPoint(uBin,Momentum,MiniumValue);
                if(CurrentValue<MiniumValue) gDownLimit.SetPoint(uBin,Momentum,CurrentValue);
                //printf("k=%.3f; Down=%.3f\n",CurrentValue);
                gUpLimit.GetPoint(uBin,Momentum,MaximumValue);
                if(CurrentValue>MaximumValue) gUpLimit.SetPoint(uBin,Momentum,CurrentValue);
                //printf("k=%.3f; Up=%.3f\n",CurrentValue);
            }

        }
        delete SystFile;
    }

//return;
    TGraphErrors *grFemtoLL;
    grFemtoLL = FemtoModelFitBandsSimple(&gDownLimit, &gUpLimit);
    grFemtoLL->SetFillColor(fColors[5]);
    grFemtoLL->SetLineColor(fColors[5]);
    grFemtoLL->SetLineWidth(5);

    TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
    Can_CF_LL->SetRightMargin(right);
    Can_CF_LL->SetTopMargin(top);

    TF1 *baselineLL = new TF1("baselineLL", "pol1", 0, 1000);
    baselineLL->SetParameter(0,parA);
    baselineLL->SetParameter(1,parB);
    baselineLL->SetLineStyle(2);
    baselineLL->SetLineColor(fFillColors[6]);
    baselineLL->SetLineWidth(4);

    TGraph *grFakeLL = new TGraph();
    grFakeLL->SetLineColor(fColors[5]);
    grFakeLL->SetLineWidth(4);

    TF1* qsLL = new TF1("qsLL","([0]+[1]*x)*(1.-[2]*0.5*exp(-4.*x*x*[3]*[3]*5.0677*5.0677*1e-6))", 0, 1000);
    qsLL->SetParameter(0, parA);
    qsLL->SetParameter(1, parB);
    if(WhichDataSet==0){
        qsLL->SetParameter(2, 0.338);//THE LAMBDA PARAMETER!!!
        qsLL->SetParameter(3, 1.152);
    }
    else if(WhichDataSet==1){
        qsLL->SetParameter(2, 0.239101);//THE LAMBDA PARAMETER!!!
        qsLL->SetParameter(3, 1.401);
    }
    else{
        qsLL->SetParameter(2, 0.2994);//THE LAMBDA PARAMETER!!!
        qsLL->SetParameter(3, 1.125);
    }


    qsLL->SetLineColor(kGreen+4);
    qsLL->SetLineStyle(5);
    qsLL->SetLineWidth(4);

    unsigned NumMomBins = DataHisto->GetNbinsX();
    double MinMom = DataHisto->GetBinLowEdge(1);
    double MaxMom = DataHisto->GetXaxis()->GetBinUpEdge(NumMomBins);
    TH1F* hGeV = new TH1F("hGeV","hGeV",NumMomBins,MinMom*0.001,MaxMom*0.001);
    for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
        hGeV->SetBinContent(uBin,DataHisto->GetBinContent(uBin));
        hGeV->SetBinError(uBin,DataHisto->GetBinError(uBin));
    }
    hGeV->SetTitle("; #it{k*} (GeV/#it{c}); #it{C}(#it{k*})");
    hGeV->GetXaxis()->SetRangeUser(0, 0.46);
    hGeV->GetXaxis()->SetNdivisions(505);
    hGeV->GetYaxis()->SetRangeUser(0.35, 2.);
    hGeV->SetFillColor(fFillColors[0]);
    SetStyleHisto(hGeV,2,0);
    hGeV->SetMarkerSize(1.5);
    //hGeV->Draw();

    DataHisto->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    DataHisto->GetXaxis()->SetRangeUser(0, 460);
    DataHisto->GetXaxis()->SetNdivisions(505);
    DataHisto->GetYaxis()->SetRangeUser(0.35, 2.);
    DataHisto->SetFillColor(fFillColors[0]);
    SetStyleHisto(DataHisto,2,0);
    DataHisto->Draw();

    TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError(DataHisto, DataSystHisto, 5);
    Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);

    //Tgraph_syserror_LL_ALAL->Draw("ap");
    baselineLL->Draw("same");
    //Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
    //Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0, 0.2);
    //Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
    //Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.35, 2.);
    if(grFemtoLL) grFemtoLL->Draw("l3 same");
    qsLL->Draw("same");
    //hGeV->Draw("same");
    DataHisto->Draw("same");

    Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror_LL_ALAL->Draw("2 same");
    //DataHisto->Draw("pe same");

    unsigned NumRows=4;
    TLegend *legLL2 = new TLegend(0.49,0.87-0.07*NumRows,0.70,0.87);//lbrt
    legLL2->SetBorderSize(0);
    legLL2->SetTextFont(42);
    legLL2->SetTextSize(gStyle->GetTextSize()*0.75);
    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hGeV->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    legLL2->AddEntry(hCk_Fake, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs", "fpe");
//  legLL2->AddEntry(hist_CF_LL_ALAL_exp[2], "with Syst. uncertainties", "");
    legLL2->AddEntry(baselineLL,"Baseline","l");
    legLL2->AddEntry(grFakeLL,"Femtoscopic fit","l");
    legLL2->AddEntry(qsLL,"Quantum statistics","l");
    legLL2->Draw("same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
    BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
    //if(WhichDataSet==0) BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 13 TeV"));
    //else if(WhichDataSet==1) BeamText.DrawLatex(0.55, 0.825, Form("p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    //else BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 7 TeV"));
    if(WhichDataSet==0) BeamText.DrawLatex(0.50, 0.90, "ALICE pp #sqrt{#it{s}} = 13 TeV");
    else if(WhichDataSet==1) BeamText.DrawLatex(0.50, 0.90, "ALICE p#minusPb #sqrt{#it{s}} = 5.02 TeV");
    else BeamText.DrawLatex(0.50, 0.90, "ALICE pp #sqrt{#it{s}} = 7 TeV");

    Can_CF_LL->Print(FitFolder+TString::Format("PlotLamLamFit_%u.pdf",WhichDataSet));

    //grFemtoLL->Draw("l3 same");

    delete grFemtoLL;
    delete Can_CF_LL;
    delete baselineLL;
    delete legLL2;
    delete grFakeLL;
    delete hGeV;
    delete hCk_Fake;
    delete qsLL;
    delete DataFile;
    delete DataSystFile;
}

void plotManyLambdaLambdaModels(const TString OutputFolder, const TString system)
{

    //const char expfile[] = "~/Results/LHC17p_fast/AnalysisResults.root";
    const TString expfile = system=="pp"?
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root":
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root";

    const TString MoritaFile = system=="pp"?
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/ComputeMorita_170119/ComputeMoritaPotentials_ALICE_pp_13TeV.root":
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/ComputeMorita_170119/ComputeMoritaPotentials_ALICE_pPb_5TeV.root";

    //int rebin = 1;
    //if(system!="pp") rebin=5;

    gStyle->SetCanvasPreferGL(1);
    const float right = 0.035;
    const float top = 0.035;

    const double energy = system=="pp"?13:5.02; // TeV

    SetStyle();

    TFile* InputFile = new TFile(expfile,"read");

    //hist_CF_LL_ALAL_exp[2]
    TH1F* hCk_LL = (TH1F*)InputFile->Get("hCkTotNormWeight");
    SetStyleHisto(hCk_LL, 0,0);
    hCk_LL->SetMarkerSize(1.5);

    TString HistoName;
    if(system=="pp") HistoName = "hCkTotNormWeight";
    else HistoName = "hCk_ReweightedMeV_1";

    TString input_sys_LL = "/home/dmihaylov/CernBox/Femto_pp13/data_old/C2totalsysLL.root";
    TFile* file_sys_LL = new TFile(input_sys_LL,"read");
    TH1F* hist_sys_LL = (TH1F*)file_sys_LL->Get("C2totalsysLL");
    hist_sys_LL->SetLineWidth(2.0);

    TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError(hCk_LL, hist_sys_LL, 5);

    TGraph *grFakeSyst = new TGraph();
    grFakeSyst->SetFillColor(fFillColors[0]);
    //grFakeSyst->SetLineColor(kBlack);
    grFakeSyst->SetLineColor(grFakeSyst->GetFillColor());
    grFakeSyst->SetMarkerStyle(kOpenCircle);
    grFakeSyst->SetMarkerSize(1.5);

    TFile* morita = new TFile(MoritaFile,"read");
    auto* grnd46 = (TGraphErrors*)morita->Get("gCkTheory_ND46");
    grnd46->SetFillColor(fColors[2]);
    grnd46->SetLineColor(fColors[2]);
    grnd46->SetLineWidth(4);
    grnd46->SetLineStyle(3);
    //auto* grnsc97f = convertInGev((TGraphErrors*)morita->Get("gCkTheory_NSC97f"));
    //grnsc97f->SetFillColor(fColors[5]);
    //grnsc97f->SetLineColor(fColors[5]);
    //grnsc97f->SetLineWidth(3);
    //grnsc97f->SetLineStyle(2);
    auto* grnf44 = (TGraphErrors*)morita->Get("gCkTheory_NF44");
    grnf44->SetFillColor(fColors[3]);
    grnf44->SetLineColor(fColors[3]);
    grnf44->SetLineWidth(4);
    grnf44->SetLineStyle(2);

    //this is not even a spline (as intended), but a weird trick to make ROOT make a smooth line
    //essentially we change the value of the first bin (that anyway lies outside of the plot), just so we get a smooth curve in the plot itself
    TGraphErrors* grnf44_spline = new TGraphErrors();
    grnf44_spline->SetName("grnf44_spline");
    unsigned NumSplinePts=10;
    grnf44_spline->Set(NumSplinePts);
    double kMin=0;
    double kMax=200;
    double kStep = (kMax-kMin)/double(NumSplinePts);
    for(unsigned uBin=0; uBin<NumSplinePts; uBin++){
        /*
        if(uBin==1){
            grnf44_spline->SetPoint(uBin,kMin+0.5*kStep+double(uBin)*kStep,1.8);
        }
        else if(uBin==2){
            grnf44_spline->SetPoint(uBin,kMin+0.5*kStep+double(uBin)*kStep,grnf44->Eval(kMin+0.5*kStep+double(uBin)*kStep)*1.0);
        }
        else grnf44_spline->SetPoint(uBin,kMin+0.5*kStep+double(uBin)*kStep,grnf44->Eval(kMin+0.5*kStep+double(uBin)*kStep));
        */

        //if(uBin==0) continue;
        if(uBin==0 && system!="pp") grnf44_spline->SetPoint(uBin,kMin+0.5*kStep+double(uBin)*kStep,2.7);
        else grnf44_spline->SetPoint(uBin,kMin+0.5*kStep+double(uBin)*kStep,grnf44->Eval(kMin+0.5*kStep+double(uBin)*kStep));
    }
    grnf44_spline->SetFillColor(fColors[3]);
    grnf44_spline->SetLineColor(fColors[3]);
    grnf44_spline->SetLineWidth(4);
    grnf44_spline->SetLineStyle(2);

    auto* grESC08 = (TGraphErrors*)morita->Get("gCkTheory_ESC08");
    grESC08->SetFillColor(fColors[4]);
    grESC08->SetLineColor(fColors[4]);
    grESC08->SetLineWidth(4);
    auto* grEhime = (TGraphErrors*)morita->Get("gCkTheory_Ehime");
    grEhime->SetFillColor(fColors[2]);
    grEhime->SetLineColor(fColors[2]);
    grEhime->SetLineWidth(4);
    auto* grHKMYY = (TGraphErrors*)morita->Get("gCkTheory_HKMYY");
    grHKMYY->SetFillColor(fColors[5]);
    grHKMYY->SetLineColor(fColors[5]);
    grHKMYY->SetLineWidth(4);
    auto* grNoSI = (TGraphErrors*)morita->Get("gCkTheory_NoSI");
    grNoSI->SetFillColor(kGreen+4);
    grNoSI->SetLineColor(kGreen+4);
    grNoSI->SetLineWidth(4);
    grNoSI->SetLineStyle(5);

    gStyle->SetErrorX(0.001);

    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.85);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.495, 0.875, "ALICE Preliminary");
    if(system=="pp") BeamText.DrawLatex(0.495, 0.825, Form("%s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
    else BeamText.DrawLatex(0.495, 0.825, Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", "p#minusPb", energy));


    TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
    Can_CF_LL->cd(0);
    Can_CF_LL->SetRightMargin(right);
    Can_CF_LL->SetTopMargin(top);
    //Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);
    //Tgraph_syserror_LL_ALAL->Draw("ap");
    //Tgraph_syserror_LL_ALAL->SetTitle("; k* (GeV/#it{c}); #it{C}(k*)");
    //Tgraph_syserror_LL_ALAL->GetXaxis()->SetRangeUser(0.35, 2.);
    //Tgraph_syserror_LL_ALAL->GetXaxis()->SetNdivisions(505);
    //Tgraph_syserror_LL_ALAL->GetYaxis()->SetRangeUser(0.35, 2.);

    hCk_LL->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hCk_LL->GetXaxis()->SetRangeUser(0, 200);
    hCk_LL->GetXaxis()->SetNdivisions(505);
    hCk_LL->GetYaxis()->SetRangeUser(0.35, 2.);
    hCk_LL->SetFillColor(fFillColors[0]);
    SetStyleHisto(hCk_LL,2,0);
    hCk_LL->Draw();


    //grnf44->Draw("C same");
    //grnf44->Draw("l same");
    //grnf44_spline->SetLineColor(kRed);
    grnf44_spline->Draw("C same");
    grEhime->Draw("Ce3 same");
    grESC08->Draw("Ce3 same");
    grnd46->Draw("C same");
    grHKMYY->Draw("Ce3 same");
    grNoSI->Draw("C same");
    //grnsc97f->Draw("l same");
    Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror_LL_ALAL->Draw("2 same");
    hCk_LL->Draw("pe same");
    unsigned NumRows2 = 1;
    TLegend *legLL2 = new TLegend(0.49,0.87-0.07*NumRows2,0.70,0.87);//lbrt
    legLL2->SetBorderSize(0);
    legLL2->SetTextFont(42);
    legLL2->SetTextSize(gStyle->GetTextSize()*0.75);
    legLL2->AddEntry(grFakeSyst, "#Lambda#minus#Lambda #oplus #bar{#Lambda}#minus#bar{#Lambda} pairs", "pef");
    legLL2->Draw("same");

    unsigned NumRows = 3;
    TLegend *legLL = new TLegend(0.49,0.80-0.07*NumRows,0.95,0.80);
    legLL->SetBorderSize(0);
    legLL->SetTextFont(42);
    legLL->SetTextSize(gStyle->GetTextSize()*0.75);
    legLL->SetNColumns(2);
    legLL->AddEntry(grnd46, "ND46", "l");
    //legLL->AddEntry(grnsc97f, "NSC97f", "l");
    legLL->AddEntry(grnf44, "NF44", "l");
    legLL->AddEntry(grEhime, "Ehime", "l");
    legLL->AddEntry(grESC08, "ESC08", "l");
    legLL->AddEntry(grHKMYY, "HKMYY", "l");
    legLL->AddEntry(grNoSI, "Quantum statistics", "l");
    legLL->Draw("same");
    //ref.DrawLatex(0.555, 0.63, "PRC 91 (2015) 024916.");
    if(system=="pp") BeamText.DrawLatex(0.50, 0.90, Form("ALICE %s #sqrt{#it{s}} = %.0f TeV", system.Data(), energy));
    else BeamText.DrawLatex(0.50, 0.90, Form("ALICE %s #sqrt{#it{s}_{NN}} = %.2f TeV", "p#minusPb", energy));
    TDatime DateTime;
    //Can_CF_LL->Print(OutputFolder+TString::Format("%i-%i-%i-CF_LL_Models_%s.pdf",DateTime.GetYear(),DateTime.GetMonth(),DateTime.GetDay(),system.Data()));
    Can_CF_LL->Print(OutputFolder+TString::Format("CF_LL_Models_%s.pdf",system.Data()));

    delete file_sys_LL;
    delete grnf44_spline;
    delete InputFile;
}

void Friendship(const TString& FolderName, const TString& InputFileName) {
  // Data
  //TFile *file2 = TFile::Open("CombinedMap_ForIlya.root");
  TFile *file2 = TFile::Open(FolderName+InputFileName);
  TH2F *sigma, *ledniSucks;
  ledniSucks = (TH2F*)file2->Get("hGlobMinCk_Def");
  int nBinsX=ledniSucks->GetXaxis()->GetNbins();
  float xMin=ledniSucks->GetXaxis()->GetBinLowEdge(1);
  float xMax=ledniSucks->GetXaxis()->GetBinUpEdge(nBinsX);
  int nBinsY=ledniSucks->GetYaxis()->GetNbins();
  float yMin=ledniSucks->GetYaxis()->GetBinLowEdge(1);
  float yMax=ledniSucks->GetYaxis()->GetBinUpEdge(nBinsY);
  std::cout << nBinsX << '\t' << xMin << '\t' << xMax << '\t' <<
      nBinsY << '\t' << yMin << '\t' << yMax << std::endl;
  TH2F* invLedniSucks=new TH2F("invLedniSucks","invLedniSucks",nBinsX,xMin,xMax
                               ,nBinsY,yMin,yMax);
  for (int iXbin=0;iXbin<=nBinsX;++iXbin) {
    for (int iYbin=0;iYbin<=nBinsY;++iYbin) {
      if (ledniSucks->GetBinContent(iXbin,iYbin)==0) {
        invLedniSucks->SetBinContent(iXbin,iYbin,15);
      } else {
        invLedniSucks->SetBinContent(iXbin,iYbin,0);
      }
    }
  }

  TFile *output=new TFile(FolderName+"Friendship.root","RECREATE");
  invLedniSucks->Write();
}


//from DataSample,SourceType etc. I can get the fit result file, which also has in it the data
//=> the data file name is this file.
//the graph with the fit result in the same file is to be used for default, up and lower limit!
void Plot_pL_FAST(      const TString& WorkFolder,
                        const TString& DataSample, const TString& SourceType, const TString& FittingMode_pp, const TString& FittingMode_pL,
                        const TString& pL_Model1, const TString& pL_Model2,
                        const TString& LegendSource_line1, const TString& LegendSource_line2,
                        const double parA, const double parB, const TString DataSystFileName){

                //const TString FitFolder, const unsigned WhichDataSet, const unsigned NumSystIter,
                //   const TString DataFileName, const TString DataHistoName, const double parA, const double parB,
                //   const TString DataSystFileName){


    TString DataFileName = WorkFolder+TString::Format("Output_%s_%s_%s.root",DataSample.Data(),SourceType.Data(),FittingMode_pp.Data());
    TString DataHistoName = "hCk_ReweightedMeV_2";
    TString FitResultHistoName1 = "pL_AV18_"+pL_Model1+"_pXim_HALQCD1";
    TString FitResultHistoName2 = "pL_AV18_"+pL_Model2+"_pXim_HALQCD1";

    TFile* DataFile = new TFile(DataFileName,"read");
    if(!DataFile){
        printf("Problem with %s\n",DataFileName.Data());
        return;
    }
    TH1F* DataHisto = (TH1F*)DataFile->Get(DataHistoName);
    if(!DataHisto){
        printf("Problem with %s\n",DataHistoName.Data());
        return;
    }
    TGraph* FitResult1 = (TGraph*)DataFile->Get(FitResultHistoName1);
    TGraph* FitResult2 = (TGraph*)DataFile->Get(FitResultHistoName2);

    TGraph gDownLimit;
    TGraph gUpLimit;

    TFile* DataSystFile = new TFile(DataSystFileName);
    TH1F* DataSystHisto = (TH1F*)DataSystFile->Get("SystErrRel");
    TF1* ExtrapolOLD = NULL;
    if(!DataSystHisto){
        ExtrapolOLD = (TF1*)DataSystFile->Get("RelSysPLUnbinned");
    }
//printf("DataSystFile=%p\n",DataSystFile);
//printf("DataSystHisto=%p\n",DataSystHisto);
//printf("ExtrapolOLD=%p\n",ExtrapolOLD);

    gStyle->SetCanvasPreferGL(1);
    SetStyle();
    const float right = 0.025;
    const float top = 0.025;

    TGraphErrors *grFemtoLL_1;
    grFemtoLL_1 = FemtoModelFitBandsSimple(FitResult1, FitResult1);
    grFemtoLL_1->SetFillColor(fColors[3]);
    grFemtoLL_1->SetLineColor(fColors[3]);
    grFemtoLL_1->SetLineWidth(5);

    TGraphErrors *grFemtoLL_2;
    grFemtoLL_2 = FemtoModelFitBandsSimple(FitResult2, FitResult2);
    grFemtoLL_2->SetFillColor(fColors[1]);
    grFemtoLL_2->SetLineColor(fColors[1]);
    grFemtoLL_2->SetLineWidth(5);

    TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
    Can_CF_LL->SetRightMargin(right);
    Can_CF_LL->SetTopMargin(top);

    TF1 *baselineLL = new TF1("baselineLL", "pol1", 0, 1000);
    baselineLL->SetParameter(0,parA);
    baselineLL->SetParameter(1,parB);
    baselineLL->SetLineStyle(2);
    baselineLL->SetLineColor(fFillColors[6]);
    baselineLL->SetLineWidth(4);

    TGraph *grFakeLL = new TGraph();
    grFakeLL->SetLineColor(fColors[5]);
    grFakeLL->SetLineWidth(4);

    unsigned NumMomBins = DataHisto->GetNbinsX();
    double MinMom = DataHisto->GetBinLowEdge(1);
    double MaxMom = DataHisto->GetXaxis()->GetBinUpEdge(NumMomBins);
    TH1F* hGeV = new TH1F("hGeV","hGeV",NumMomBins,MinMom*0.001,MaxMom*0.001);
    for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
        hGeV->SetBinContent(uBin,DataHisto->GetBinContent(uBin));
        hGeV->SetBinError(uBin,DataHisto->GetBinError(uBin));
    }
    hGeV->SetTitle("; #it{k*} (GeV/#it{c}); #it{C}(#it{k*})");
    hGeV->GetXaxis()->SetRangeUser(0, 0.32);
    hGeV->GetXaxis()->SetNdivisions(505);
    hGeV->GetYaxis()->SetRangeUser(0.85, 2.3);
    hGeV->SetFillColor(fFillColors[0]);
    SetStyleHisto(hGeV,2,0);
    hGeV->SetMarkerSize(1.5);
    //hGeV->Draw();

    DataHisto->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    DataHisto->GetXaxis()->SetRangeUser(0, 320);
    DataHisto->GetXaxis()->SetNdivisions(505);
    DataHisto->GetYaxis()->SetRangeUser(0.85, 2.3);
    DataHisto->SetFillColor(fFillColors[0]);
    SetStyleHisto(DataHisto,2,0);
    DataHisto->Draw();

    TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError_FAST(DataHisto, DataSystHisto, ExtrapolOLD, 3);
    Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);

    //baselineLL->Draw("same");

    if(grFemtoLL_1) grFemtoLL_1->Draw("l3 same");
    if(grFemtoLL_2) grFemtoLL_2->Draw("l3 same");
    //hGeV->Draw("same");
    DataHisto->Draw("same");

    Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror_LL_ALAL->Draw("2 same");
    //DataHisto->Draw("pe same");

    unsigned NumRows=4;
    TLegend *legLL2 = new TLegend(0.49,0.76-0.04*NumRows,0.73,0.76);//lbrt
    legLL2->SetBorderSize(0);
    legLL2->SetTextFont(42);
    legLL2->SetTextSize(gStyle->GetTextSize()*0.70);
    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hGeV->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    legLL2->AddEntry(hCk_Fake, "p#minus#Lambda #oplus #bar{p}#minus#bar{#Lambda} pairs", "fpe");
//  legLL2->AddEntry(hist_CF_LL_ALAL_exp[2], "with Syst. uncertainties", "");
    //legLL2->AddEntry(baselineLL,"Baseline","l");
    legLL2->AddEntry(grFemtoLL_1,"Femtoscopic fit (#chiEFT LO)","l");
    legLL2->AddEntry(grFemtoLL_2,"Femtoscopic fit (#chiEFT NLO)","l");
    legLL2->Draw("same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.80);
    BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
    //if(WhichDataSet==0) BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 13 TeV"));
    //else if(WhichDataSet==1) BeamText.DrawLatex(0.55, 0.825, Form("p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    //else BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 7 TeV"));
    BeamText.DrawLatex(0.50, 0.91, "ALICE Preliminary");
    if(DataSample=="pp13TeV_MB_Run2paper") BeamText.DrawLatex(0.50, 0.86, "pp #sqrt{#it{s}} = 13 TeV");
    else if(DataSample.Contains("pPb")&&DataSample.Contains("5TeV")) BeamText.DrawLatex(0.50, 0.86, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    else if(DataSample=="pp13TeV_HM_March19") BeamText.DrawLatex(0.50, 0.86, "pp (HM) #sqrt{#it{s}} = 13 TeV");
    else BeamText.DrawLatex(0.50, 0.86, "ALICE pp #sqrt{#it{s}} = 7 TeV");

    TLatex BeamTextSource;
    BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.70);
    BeamTextSource.SetNDC(kTRUE);
    BeamTextSource.DrawLatex(0.50, 0.81, LegendSource_line1);



//INLET -------------------------------------------------------------------------------------------------------------------

    TH1F* DataHisto_Inlet = (TH1F*)DataHisto->Clone("DataHisto_Inlet");
    DataHisto_Inlet->SetMarkerSize(DataHisto->GetMarkerSize()*0.67);
    DataHisto_Inlet->SetLineWidth(DataHisto->GetLineWidth()*0.67);
    DataHisto_Inlet->GetXaxis()->SetTitleSize(DataHisto->GetXaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetLabelSize(DataHisto->GetXaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetRangeUser(120, 320);
    DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

    DataHisto_Inlet->GetYaxis()->SetTitleSize(DataHisto->GetYaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetLabelSize(DataHisto->GetYaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetTitleOffset(DataHisto->GetYaxis()->GetTitleOffset()*0.67);
    DataHisto_Inlet->GetYaxis()->SetRangeUser(0.98, 1.045);

    TGraph* grFemtoLL_1_Inlet = (TGraph*)grFemtoLL_1->Clone("grFemtoLL_1_Inlet");
    grFemtoLL_1_Inlet->SetLineWidth(grFemtoLL_1->GetLineWidth()*0.67);

    TGraph* grFemtoLL_2_Inlet = (TGraph*)grFemtoLL_2->Clone("grFemtoLL_2_Inlet");
    grFemtoLL_2_Inlet->SetLineWidth(grFemtoLL_2->GetLineWidth()*0.67);

    const double fXMinInlet=0.35;
    const double fYMinInlet=0.25;
    const double fXMaxInlet=0.95;
    const double fYMaxInlet=0.57;
    TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                             fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    DataHisto_Inlet->Draw();
    if(grFemtoLL_1) grFemtoLL_1_Inlet->Draw("l3 same");
    if(grFemtoLL_2) grFemtoLL_2_Inlet->Draw("l3 same");
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror_LL_ALAL->Draw("2 same");

    Can_CF_LL->Print(WorkFolder+TString::Format("PlotLamLamFit_%s_%s_%s.pdf",DataSample.Data(),SourceType.Data(),FittingMode_pp.Data()));

    delete inset_pad;
    delete DataHisto_Inlet;
    delete grFemtoLL_1_Inlet;
    delete grFemtoLL_2_Inlet;

    delete grFemtoLL_1;
    delete grFemtoLL_2;
    delete Can_CF_LL;
    delete baselineLL;
    delete legLL2;
    delete grFakeLL;
    delete hGeV;
    delete hCk_Fake;

    delete DataFile;
    delete DataSystFile;
}


//from DataSample,SourceType etc. I can get the fit result file, which also has in it the data
//=> the data file name is this file.
//the graph with the fit result in the same file is to be used for default, up and lower limit!


void Plot_pL_FASTsyst( const int& WhichBaseline

                      //const TString& DataFileName, const TString& DataHistoName,
                      //  const TString&

                      //const TString& WorkFolder,
                        //const TString& DataSample, const TString& SourceType, const TString& FittingMode_pp, const TString& FittingMode_pL,
                        //const TString& pL_Model1, const TString& pL_Model2,
                        //const TString& LegendSource_line1, const TString& LegendSource_line2,
                        //const double parA, const double parB, const TString DataSystFileName
                        ){

                //const TString FitFolder, const unsigned WhichDataSet, const unsigned NumSystIter,
                //   const TString DataFileName, const TString DataHistoName, const double parA, const double parB,
                //   const TString DataSystFileName){


//MODE = 0 -> the coupled fit
//MODE = 1 -> the one without coupled
    int MODE = 0;
    double MinMomentum=0;
    double MaxMomentum=336;
    //int DefaultIterLO = 32400;
    //int DefaultIterNLO = 93960;
    int DefaultIterLO = 346274;
    int DefaultIterNLO = 415529;

    //we take 10HM not to redo further plots
    TString DataFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/CFOutput_pL.root";
    //TString DataFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample12HM/CFOutput_pL_0.root";
    TString DataHistoName = "hCk_ReweightedMeV_2";
    TFile* DataFile = new TFile(DataFileName,"read");
    TH1F* hData = (TH1F*)DataFile->Get(DataHistoName);

    TString DataSystFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample12HM/Systematics_pL.root";
    TString DataSystHistoName = "hRelSyst";//"RelSysPLUnbinned" "hRelSyst" "SystErrRel"
    TFile* SystFile = new TFile(DataSystFileName,"read");
    TH1F* hSyst = (TH1F*)SystFile->Get(DataSystHistoName);
    TF1* fSyst = (TF1*)SystFile->Get(DataSystHistoName);
//SystematicsAdd_100419_2 pm15
//SystematicsAdd_120419 pm10
//SystematicsAddMeson_230419
//SystematicsAdd_310519
//SystematicsAdd_010519
    TString FitSystFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_090519/";
    TString FitSystFileName = FitSystFolder+"NTfile_0.root";

    TFile* ntFile = new TFile(FitSystFileName,"read");
    TNtuple* ntResult = (TNtuple*)ntFile->Get("ntResult");
    unsigned NumNtEntries = ntResult->GetEntries();

    Float_t Iter;
    Float_t Config;
    Float_t Data;
    Float_t SourceType;
    Float_t SourceScale;
    Float_t Potential;
    Float_t Baseline;
    Float_t FemRan;
    Float_t FitRan;
    Float_t pFrac;
    Float_t LamFrac;
    Float_t kcVar;
    Float_t mTbin;
    Float_t FemtoMin;
    Float_t FemtoMax;
    Float_t BlMin;
    Float_t BlMax;
    Float_t p_a;
    Float_t e_a;
    Float_t p_b;
    Float_t e_b;
    Float_t p_c;
    Float_t e_c;
    Float_t p_Cl;
    Float_t e_Cl;
    Float_t p_kc;
    Float_t e_kc;
    Float_t p_sor0;
    Float_t e_sor0;
    Float_t p_sor1;
    Float_t e_sor1;
    Float_t chi2;
    Float_t ndf;
    Float_t chi2_312;
    Float_t ndf_312;

    ntResult->SetBranchAddress("Iter",&Iter);
    ntResult->SetBranchAddress("Config",&Config);
    ntResult->SetBranchAddress("Data",&Data);
    ntResult->SetBranchAddress("SourceType",&SourceType);
    ntResult->SetBranchAddress("SourceScale",&SourceScale);
    ntResult->SetBranchAddress("Potential",&Potential);
    ntResult->SetBranchAddress("Baseline",&Baseline);
    ntResult->SetBranchAddress("FemRan",&FemRan);
    ntResult->SetBranchAddress("FitRan",&FitRan);
    ntResult->SetBranchAddress("pFrac",&pFrac);
    ntResult->SetBranchAddress("LamFrac",&LamFrac);
    ntResult->SetBranchAddress("kcVar",&kcVar);
    ntResult->SetBranchAddress("mTbin",&mTbin);
    ntResult->SetBranchAddress("FemtoMin",&FemtoMin);
    ntResult->SetBranchAddress("FemtoMax",&FemtoMax);
    ntResult->SetBranchAddress("BlMin",&BlMin);
    ntResult->SetBranchAddress("BlMax",&BlMax);
    ntResult->SetBranchAddress("p_a",&p_a);
    ntResult->SetBranchAddress("e_a",&e_a);
    ntResult->SetBranchAddress("p_b",&p_b);
    ntResult->SetBranchAddress("e_b",&e_b);
    ntResult->SetBranchAddress("p_c",&p_c);
    ntResult->SetBranchAddress("e_c",&e_c);
    ntResult->SetBranchAddress("p_Cl",&p_Cl);
    ntResult->SetBranchAddress("e_Cl",&e_Cl);
    ntResult->SetBranchAddress("p_kc",&p_kc);
    ntResult->SetBranchAddress("e_kc",&e_kc);
    ntResult->SetBranchAddress("p_sor0",&p_sor0);
    ntResult->SetBranchAddress("e_sor0",&e_sor0);
    ntResult->SetBranchAddress("p_sor1",&p_sor1);
    ntResult->SetBranchAddress("e_sor1",&e_sor1);
    ntResult->SetBranchAddress("chi2",&chi2);
    ntResult->SetBranchAddress("ndf",&ndf);
    ntResult->SetBranchAddress("chi2_312",&chi2_312);
    ntResult->SetBranchAddress("ndf_312",&ndf_312);

    TGraph gLowerLO;
    gLowerLO.SetName("gLowerLO");
    TGraph gUpperLO;
    gUpperLO.SetName("gUpperLO");
    TGraph gBestLO;
    gBestLO.SetName("gBestLO");

    //baseline
    TGraph bOuterLowerLO;
    bOuterLowerLO.SetName("bOuterLowerLO");
    TGraph bOuterUpperLO;
    bOuterUpperLO.SetName("bOuterUpperLO");
    TGraph bOuterBestLO;
    bOuterBestLO.SetName("bOuterBestLO");
    TGraph bLowerLO;
    bLowerLO.SetName("bLowerLO");
    TGraph bUpperLO;
    bUpperLO.SetName("bUpperLO");
    TGraph bBestLO;
    bBestLO.SetName("bBestLO");

    TGraph gLowerNLO;
    gLowerNLO.SetName("gLowerNLO");
    TGraph gUpperNLO;
    gUpperNLO.SetName("gUpperNLO");
    TGraph gBestNLO;
    gBestNLO.SetName("gBestNLO");

    TGraph bOuterLowerNLO;
    bOuterLowerNLO.SetName("bOuterLowerNLO");
    TGraph bOuterUpperNLO;
    bOuterUpperNLO.SetName("bOuterUpperNLO");
    TGraph bOuterBestNLO;
    bOuterBestNLO.SetName("bOuterBestNLO");
    TGraph bLowerNLO;
    bLowerNLO.SetName("bLowerNLO");
    TGraph bUpperNLO;
    bUpperNLO.SetName("bUpperNLO");
    TGraph bBestNLO;
    bBestNLO.SetName("bBestNLO");

    double BestChi2NdfLO = 1e6;
    double WorstChi2NdfLO = 0;
    double BestChi2NdfNLO = 1e6;
    double WorstChi2NdfNLO = 0;
    unsigned NumAcceptedEntries = 0;

    double GlobalBestPval_LO = 0;
    double GlobalBestPval_NLO = 0;

    TString LegendSource_line1 = "Gaussian core with m_{T} scaling";
    //TString LegendSource_line1 = "Levy core with m_{T} scaling";
    //TString LegendSource_line1 = "Gaussian source";

    TH1F* hChi2NdfLO = new TH1F("hChi2NdfLO","hChi2NdfLO",1024,0,64);
    TH1F* hChi2NdfNLO = new TH1F("hChi2NdfNLO","hChi2NdfNLO",1024,0,64);
    TH1F* hNsigmaLO = new TH1F("hNsigmaLO","hNsigmaLO",1024,0,64);
    TH1F* hNsigmaNLO = new TH1F("hNsigmaNLO","hNsigmaNLO",1024,0,64);

    //bin-by-bin deviation in nSigma
    //we use TH2F to save all variations
    //than we take the mean and stdv to define upper and lower limit
    const unsigned DevPlotNumBins = 47;
    TH2F* hVarDeviationLO = new TH2F("hVarDeviationLO","hVarDeviationLO",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,1024,-40,40);
    TH2F* hVarDeviationNLO = new TH2F("hVarDeviationNLO","hVarDeviationNLO",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,1024,-40,40);

    //TH2F* hVarDevRatioNsigmaLO = new TH2F("hVarDevRatioNsigmaLO","hVarDevRatioNsigmaLO",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,128,0.85,1.15);
    //TH2F* hVarDevRatioNsigmaNLO = new TH2F("hVarDevRatioNsigmaNLO","hVarDevRatioNsigmaNLO",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,128,0.85,1.15);

    TH2F* hVarDevRatioLO = new TH2F("hVarDevRatioLO","hVarDevRatioLO",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,128,0.85,1.15);
    TH2F* hVarDevRatioNLO = new TH2F("hVarDevRatioNLO","hVarDevRatioNLO",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,128,0.85,1.15);

    TGraphErrors gVarDeviationLO_Stdv;
    gVarDeviationLO_Stdv.SetName("gVarDeviationLO_Stdv");
    gVarDeviationLO_Stdv.SetLineWidth(0);
    gVarDeviationLO_Stdv.SetLineColor(kGreen+3);
    gVarDeviationLO_Stdv.SetMarkerSize(0);
    gVarDeviationLO_Stdv.SetMarkerStyle(20);
    gVarDeviationLO_Stdv.SetMarkerColor(kGreen+3);
    gVarDeviationLO_Stdv.SetFillColor(kGreen+3);

    TGraphErrors gVarDeviationNLO_Stdv;
    gVarDeviationNLO_Stdv.SetName("gVarDeviationNLO_Stdv");
    gVarDeviationNLO_Stdv.SetLineWidth(0);
    gVarDeviationNLO_Stdv.SetLineColor(kRed+1);
    gVarDeviationNLO_Stdv.SetMarkerSize(0);
    gVarDeviationNLO_Stdv.SetMarkerStyle(20);
    gVarDeviationNLO_Stdv.SetMarkerColor(kRed+1);
    gVarDeviationNLO_Stdv.SetFillColor(kRed+1);


    TString LegendSource_line2;
    if(WhichBaseline==0) LegendSource_line2 = "Constant baseline";
    else if(WhichBaseline==11) LegendSource_line2 = "Linear baseline";
    else if(WhichBaseline==12) LegendSource_line2 = "Quadratic baseline";

    double DefaultChi2NdfLO=0;
    double DefaultNsigmaLO=0;
    double DefaultChi2NdfNLO=0;
    double DefaultNsigmaNLO=0;

    //for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
    //    ntFile->cd();
    //    ntResult->GetEntry(uEntry);
    //}

    for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
        ntFile->cd();
        ntResult->GetEntry(uEntry);
        int ITER = int(Iter);
        int CONFIG = int(Config);
        int Id_LO;
        int Id_NLO;
        if(CONFIG<10){
            Id_LO = 1;
            Id_NLO = 11;
        }
        else if(CONFIG<100){
            Id_LO = 0;
            Id_NLO = 10;
        }
        else if(CONFIG<1000){
            Id_LO = 1012;
            Id_NLO = 1013;
        }

//printf("Iter=%i; CONFIG=%i; LO=%i; NLO=%i\n",ITER,CONFIG,Id_LO,Id_NLO);
        //! Conditions

        //if(CONFIG!=0&&CONFIG!=2) continue; //select Gauss+Reso and Gauss
        if(CONFIG!=0) continue; //select Gauss+Reso
        //if(CONFIG!=1) continue; //select Levy+Reso
        //if(CONFIG!=2) continue; //select Gauss
        if(round(Baseline)!=WhichBaseline) continue;
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        if(Id_LO==Potential){
            hChi2NdfLO->Fill(chi2_312/ndf_312);
            hNsigmaLO->Fill(sqrt(2)*TMath::ErfcInverse(TMath::Prob(chi2_312,int(ndf_312))));
            if(DefaultIterLO==int(Iter)){
                DefaultChi2NdfLO = chi2_312/ndf_312;
                DefaultNsigmaLO = sqrt(2)*TMath::ErfcInverse(TMath::Prob(chi2_312,int(ndf_312)));
            }
        }
        else if(Id_NLO==Potential){
            hChi2NdfNLO->Fill(chi2_312/ndf_312);
            hNsigmaNLO->Fill(sqrt(2)*TMath::ErfcInverse(TMath::Prob(chi2_312,int(ndf_312))));
            if(DefaultIterNLO==int(Iter)){
                DefaultChi2NdfNLO = chi2_312/ndf_312;
                DefaultNsigmaNLO = sqrt(2)*TMath::ErfcInverse(TMath::Prob(chi2_312,int(ndf_312)));
            }
        }
        else{
            printf("PROBLEM\n");
        }
    }

    double MeanChi2NdfLO = hChi2NdfLO->GetMean();
    double StdvChi2NdfLO = hChi2NdfLO->GetStdDev();
    double LenChi2NdfLO;
    double UniformStdvChi2NdfLO;
    double MinChi2NdfLO=1e6;
    double MaxChi2NdfLO=-1e6;

    double MeanChi2NdfNLO = hChi2NdfNLO->GetMean();
    double StdvChi2NdfNLO = hChi2NdfNLO->GetStdDev();
    double LenChi2NdfNLO;
    double UniformStdvChi2NdfNLO;
    double MinChi2NdfNLO=1e6;
    double MaxChi2NdfNLO=-1e6;

    double MeanNsigmaLO = hNsigmaLO->GetMean();
    double StdvNsigmaLO = hNsigmaLO->GetStdDev();
    double LenNsigmaLO;
    double UniformStdvNsigmaLO;
    double MinNsigmaLO=1e6;
    double MaxNsigmaLO=-1e6;

    double MeanNsigmaNLO = hNsigmaNLO->GetMean();
    double StdvNsigmaNLO = hNsigmaNLO->GetStdDev();
    double LenNsigmaNLO;
    double UniformStdvNsigmaNLO;
    double MinNsigmaNLO=1e6;
    double MaxNsigmaNLO=-1e6;

    for(unsigned uBin=0; uBin<hChi2NdfLO->GetNbinsX(); uBin++){
        if(hChi2NdfLO->GetBinContent(uBin+1)&&MinChi2NdfLO>hChi2NdfLO->GetBinCenter(uBin+1)){MinChi2NdfLO=hChi2NdfLO->GetBinCenter(uBin+1);}
        if(hChi2NdfLO->GetBinContent(uBin+1)&&MaxChi2NdfLO<hChi2NdfLO->GetBinCenter(uBin+1)){MaxChi2NdfLO=hChi2NdfLO->GetBinCenter(uBin+1);}
    }
    LenChi2NdfLO = MaxChi2NdfLO-MinChi2NdfLO;
    UniformStdvChi2NdfLO = LenChi2NdfLO/sqrt(12.);

    for(unsigned uBin=0; uBin<hChi2NdfNLO->GetNbinsX(); uBin++){
        if(hChi2NdfNLO->GetBinContent(uBin+1)&&MinChi2NdfNLO>hChi2NdfNLO->GetBinCenter(uBin+1)){MinChi2NdfNLO=hChi2NdfNLO->GetBinCenter(uBin+1);}
        if(hChi2NdfNLO->GetBinContent(uBin+1)&&MaxChi2NdfNLO<hChi2NdfNLO->GetBinCenter(uBin+1)){MaxChi2NdfNLO=hChi2NdfNLO->GetBinCenter(uBin+1);}
    }
    LenChi2NdfNLO = MaxChi2NdfNLO-MinChi2NdfNLO;
    UniformStdvChi2NdfNLO = LenChi2NdfNLO/sqrt(12.);

    for(unsigned uBin=0; uBin<hNsigmaLO->GetNbinsX(); uBin++){
        if(hNsigmaLO->GetBinContent(uBin+1)&&MinNsigmaLO>hNsigmaLO->GetBinCenter(uBin+1)){MinNsigmaLO=hNsigmaLO->GetBinCenter(uBin+1);}
        if(hNsigmaLO->GetBinContent(uBin+1)&&MaxNsigmaLO<hNsigmaLO->GetBinCenter(uBin+1)){MaxNsigmaLO=hNsigmaLO->GetBinCenter(uBin+1);}
    }
    LenNsigmaLO = MaxNsigmaLO-MinNsigmaLO;
    UniformStdvNsigmaLO = LenNsigmaLO/sqrt(12.);

    for(unsigned uBin=0; uBin<hNsigmaNLO->GetNbinsX(); uBin++){
        if(hNsigmaNLO->GetBinContent(uBin+1)&&MinNsigmaNLO>hNsigmaNLO->GetBinCenter(uBin+1)){MinNsigmaNLO=hNsigmaNLO->GetBinCenter(uBin+1);}
        if(hNsigmaNLO->GetBinContent(uBin+1)&&MaxNsigmaNLO<hNsigmaNLO->GetBinCenter(uBin+1)){MaxNsigmaNLO=hNsigmaNLO->GetBinCenter(uBin+1);}
    }
    LenNsigmaNLO = MaxNsigmaNLO-MinNsigmaNLO;
    UniformStdvNsigmaNLO = LenNsigmaNLO/sqrt(12.);

    double BestNsigma_LO = sqrt(2)*TMath::ErfcInverse(GlobalBestPval_LO);
    double BestNsigma_NLO = sqrt(2)*TMath::ErfcInverse(GlobalBestPval_NLO);

    double LowCi_Chi2NdfLO;
    double UpCi_Chi2NdfLO;
    double Median_Chi2NdfLO = GetCentralInterval(*hChi2NdfLO,0.68,LowCi_Chi2NdfLO,UpCi_Chi2NdfLO,true);

    double LowCi_Chi2NdfNLO;
    double UpCi_Chi2NdfNLO;
    double Median_Chi2NdfNLO = GetCentralInterval(*hChi2NdfNLO,0.68,LowCi_Chi2NdfNLO,UpCi_Chi2NdfNLO,true);

    double LowCi_NsigmaLO;
    double UpCi_NsigmaLO;
    double Median_NsigmaLO = GetCentralInterval(*hNsigmaLO,0.68,LowCi_NsigmaLO,UpCi_NsigmaLO,true);

    double LowCi_NsigmaNLO;
    double UpCi_NsigmaNLO;
    double Median_NsigmaNLO = GetCentralInterval(*hNsigmaNLO,0.68,LowCi_NsigmaNLO,UpCi_NsigmaNLO,true);

    //THANK YOU ROOT, THANKS!!!
    //so if we read the same TNtuple twice, the second iteration is SUPER slow (like c.a. x1000 times).
    //to fix this: delete the TNtuple and open it again
    delete ntResult;
    ntResult = (TNtuple*)ntFile->Get("ntResult");

    ntResult->SetBranchAddress("Iter",&Iter);
    ntResult->SetBranchAddress("Config",&Config);
    ntResult->SetBranchAddress("Data",&Data);
    ntResult->SetBranchAddress("SourceType",&SourceType);
    ntResult->SetBranchAddress("SourceScale",&SourceScale);
    ntResult->SetBranchAddress("Potential",&Potential);
    ntResult->SetBranchAddress("Baseline",&Baseline);
    ntResult->SetBranchAddress("FemRan",&FemRan);
    ntResult->SetBranchAddress("FitRan",&FitRan);
    ntResult->SetBranchAddress("pFrac",&pFrac);
    ntResult->SetBranchAddress("LamFrac",&LamFrac);
    ntResult->SetBranchAddress("kcVar",&kcVar);
    ntResult->SetBranchAddress("mTbin",&mTbin);
    ntResult->SetBranchAddress("FemtoMin",&FemtoMin);
    ntResult->SetBranchAddress("FemtoMax",&FemtoMax);
    ntResult->SetBranchAddress("BlMin",&BlMin);
    ntResult->SetBranchAddress("BlMax",&BlMax);
    ntResult->SetBranchAddress("p_a",&p_a);
    ntResult->SetBranchAddress("e_a",&e_a);
    ntResult->SetBranchAddress("p_b",&p_b);
    ntResult->SetBranchAddress("e_b",&e_b);
    ntResult->SetBranchAddress("p_c",&p_c);
    ntResult->SetBranchAddress("e_c",&e_c);
    ntResult->SetBranchAddress("p_Cl",&p_Cl);
    ntResult->SetBranchAddress("e_Cl",&e_Cl);
    ntResult->SetBranchAddress("p_kc",&p_kc);
    ntResult->SetBranchAddress("e_kc",&e_kc);
    ntResult->SetBranchAddress("p_sor0",&p_sor0);
    ntResult->SetBranchAddress("e_sor0",&e_sor0);
    ntResult->SetBranchAddress("p_sor1",&p_sor1);
    ntResult->SetBranchAddress("e_sor1",&e_sor1);
    ntResult->SetBranchAddress("chi2",&chi2);
    ntResult->SetBranchAddress("ndf",&ndf);
    ntResult->SetBranchAddress("chi2_312",&chi2_312);//chi2<=312 MeV
    ntResult->SetBranchAddress("ndf_312",&ndf_312);//number of data points up to 312 MeV

//TRandom3 rangen(11);
//printf("NumNtEntries=%u\n",NumNtEntries);
    for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
        //if(uEntry%100==0) printf("uEntry=%u\n",uEntry);
//if(rangen.Uniform()>0.01) continue;
        ntFile->cd();
        ntResult->GetEntry(uEntry);

        int ITER = int(Iter);
        int CONFIG = int(Config);
        int Id_LO;
        int Id_NLO;
        if(CONFIG<10){
            Id_LO = 1;
            Id_NLO = 11;
        }
        else if(CONFIG<100){
            Id_LO = 0;
            Id_NLO = 10;
        }
        else if(CONFIG<1000){
            Id_LO = 1012;
            Id_NLO = 1013;
        }
//printf("Iter=%i; CONFIG=%i; LO=%i; NLO=%i\n",ITER,CONFIG,Id_LO,Id_NLO);
        //! Conditions
//if(rangen.Uniform()>0.01&&ITER!=DefaultIterLO&&ITER!=DefaultIterNLO) continue;
        //if(CONFIG!=0&&CONFIG!=2) continue; //select Gauss+Reso and Gauss
        if(CONFIG!=0) continue; //select Gauss+Reso
        //if(CONFIG!=1) continue; //select Levy+Reso
        //if(CONFIG!=2) continue; //select Gauss
        if(round(Baseline)!=WhichBaseline) continue;
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        if(Potential==Id_LO&&chi2_312/double(ndf_312)<LowCi_Chi2NdfLO&&chi2_312/double(ndf_312)>UpCi_Chi2NdfLO) continue;
        if(Potential==Id_NLO&&chi2_312/double(ndf_312)<LowCi_Chi2NdfNLO&&chi2_312/double(ndf_312)>UpCi_Chi2NdfNLO) continue;
        //if(Data!=0) continue;

        bool BestSoFarLO=false;
        bool BestSoFarNLO=false;

        //we want to plot something that makes sense (i.e. is fixed to the data we plot)
        if(Data==0){
            if(Potential==Id_LO&&(chi2_312/ndf_312)>WorstChi2NdfLO&&pFrac==0&&LamFrac==0){
                WorstChi2NdfLO=chi2_312/ndf_312;
            }
            if(Potential==Id_LO&&(chi2_312/ndf_312)<BestChi2NdfLO&&pFrac==0&&LamFrac==0){
                BestChi2NdfLO=chi2_312/ndf_312;
                BestSoFarLO=true;
            }
            if(Potential==Id_LO&&TMath::Prob(chi2_312,ndf_312)>GlobalBestPval_LO){
                GlobalBestPval_LO=TMath::Prob(chi2_312,ndf_312);
    //printf("GlobalBestPval_LO=%f\n",GlobalBestPval_LO);
                //BestSoFarLO=true;
            }

            //select the best solution, demanding the default lambda parameters though
            if(Potential==Id_NLO&&(chi2_312/ndf_312)>WorstChi2NdfNLO&&pFrac==0&&LamFrac==0){
                WorstChi2NdfNLO=chi2_312/ndf_312;
            }
            if(Potential==Id_NLO&&(chi2_312/ndf_312)<BestChi2NdfNLO&&pFrac==0&&LamFrac==0){
                BestChi2NdfNLO=chi2_312/ndf_312;
                BestSoFarNLO=true;
    //printf("Hello %i -> %i\n",ITER,int(Data));
            }
            if(Potential==Id_NLO&&TMath::Prob(chi2_312,ndf_312)>GlobalBestPval_NLO){
                GlobalBestPval_NLO=TMath::Prob(chi2_312,ndf_312);
                //BestSoFarNLO=true;
            }

        }

        TFile* FileGraph = new TFile(TString::Format("%sConfig%i_Iter%i.root",FitSystFolder.Data(),CONFIG,ITER),"read");
        TH1F* hDataVar = (TH1F*)FileGraph->Get("hCk_ReweightedMeV_2");
        TGraph* FitResult = (TGraph*)FileGraph->Get("FitResult_pL");
        TF1* FitBaseline = (TF1*)FileGraph->Get("fBaseline");

        double CkMaxFit;
        double kMaxFit;
        FitResult->GetPoint(FitResult->GetN()-1,kMaxFit,CkMaxFit);

        for(unsigned uBin=0; uBin<DevPlotNumBins; uBin++){
            double kValue=hVarDeviationLO->GetXaxis()->GetBinCenter(uBin+1);
            double dataValue = hDataVar->GetBinContent(uBin+1);
            double fitValue = FitResult->Eval(kValue);
            if(kValue>kMaxFit){
                fitValue = CkMaxFit;
            }
            double dataError = hDataVar->GetBinError(uBin+1);
            double Nsigma=(fitValue-dataValue)/(dataError);
            if(Baseline==0&&kValue>480){
                printf("\nkValue=%f\n",kValue);
                printf("fitValue=%f\n",fitValue);
                printf("dataValue=%f\n",dataValue);
                printf("dataError=%f\n",dataError);
                printf("Nsigma=%f\n",Nsigma);
            }
            double Ratio=dataValue/fitValue;

            if(Potential==Id_LO){
                hVarDeviationLO->Fill(kValue,Nsigma);
                hVarDevRatioLO->Fill(kValue,Ratio);
            }
            if(Potential==Id_NLO){
                hVarDeviationNLO->Fill(kValue,Nsigma);
                hVarDevRatioNLO->Fill(kValue,Ratio);
            }
        }

        double xVal,xValLowUp;
        double yVal,yValLowUp;
        TGraph& gLower = Potential==Id_LO?gLowerLO:gLowerNLO;
        TGraph& gUpper = Potential==Id_LO?gUpperLO:gUpperNLO;

        //TGraph& gDefault = Potential==1?gDefaultLO:gDefaultNLO;
        for(int iPoint=0; iPoint<FitResult->GetN(); iPoint++){
            FitResult->GetPoint(iPoint,xVal,yVal);
            if(gLower.GetN()<iPoint+1){
                gLower.SetPoint(iPoint,xVal,yVal);
            }
            else{
                gLower.GetPoint(iPoint,yValLowUp,yValLowUp);
                if(yVal<yValLowUp){
                    gLower.SetPoint(iPoint,xVal,yVal);
                }
            }
            if(gUpper.GetN()<iPoint+1){
                gUpper.SetPoint(iPoint,xVal,yVal);
            }
            else{
                gUpper.GetPoint(iPoint,yValLowUp,yValLowUp);
                if(yVal>yValLowUp){
                    gUpper.SetPoint(iPoint,xVal,yVal);
                }
            }
            //if(ITER==DefaultIter){
            //    gDefault.SetPoint(iPoint,xVal,yVal);
            //}
            if(BestSoFarLO){
                gBestLO.SetPoint(iPoint,xVal,yVal);
            }
            if(BestSoFarNLO){
                gBestNLO.SetPoint(iPoint,xVal,yVal);
            }
        }
        int CounterBlLO=0;
        int CounterBlNLO=0;
        for(int iPoint=0; iPoint<hData->GetNbinsX(); iPoint++){
            xVal = hData->GetBinCenter(iPoint+1);
            yVal = FitBaseline->Eval(xVal);
            if(Potential==Id_LO&&BlMin<=xVal&&xVal<=BlMax){
                if(bOuterLowerLO.GetN()<CounterBlLO+1){
                    bOuterLowerLO.SetPoint(CounterBlLO,xVal,yVal);
                }
                else{
                    bOuterLowerLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal<yValLowUp){
                        bOuterLowerLO.SetPoint(CounterBlLO,xVal,yVal);
                    }
                }
                if(bOuterUpperLO.GetN()<CounterBlLO+1){
                    bOuterUpperLO.SetPoint(CounterBlLO,xVal,yVal);
                }
                else{
                    bOuterUpperLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal>yValLowUp){
                        bOuterUpperLO.SetPoint(CounterBlLO,xVal,yVal);
                    }
                }
                if(BestSoFarLO){
                    bOuterBestLO.SetPoint(CounterBlLO,xVal,yVal);
                }
                CounterBlLO++;
            }
            if(Potential==Id_NLO&&BlMin<=xVal&&xVal<=BlMax){
                if(bOuterLowerNLO.GetN()<CounterBlNLO+1){
                    bOuterLowerNLO.SetPoint(CounterBlNLO,xVal,yVal);
                }
                else{
                    bOuterLowerNLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal<yValLowUp){
                        bOuterLowerNLO.SetPoint(CounterBlNLO,xVal,yVal);
                    }
                }
                if(bOuterUpperNLO.GetN()<CounterBlNLO+1){
                    bOuterUpperNLO.SetPoint(CounterBlNLO,xVal,yVal);
                }
                else{
                    bOuterUpperNLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal>yValLowUp){
                        bOuterUpperNLO.SetPoint(CounterBlNLO,xVal,yVal);
                    }
                }
                if(BestSoFarNLO){
                    bOuterBestNLO.SetPoint(CounterBlNLO,xVal,yVal);
                }
                CounterBlNLO++;
            }
        }

        CounterBlLO=0;
        CounterBlNLO=0;
        for(int iPoint=0; iPoint<hData->GetNbinsX(); iPoint++){
            xVal = hData->GetBinCenter(iPoint+1);
            yVal = FitBaseline->Eval(xVal);
            if(Potential==Id_LO){
                if(bLowerLO.GetN()<CounterBlLO+1){
                    bLowerLO.SetPoint(CounterBlLO,xVal,yVal);
                }
                else{
                    bLowerLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal<yValLowUp){
                        bLowerLO.SetPoint(CounterBlLO,xVal,yVal);
                    }
                }
                if(bUpperLO.GetN()<CounterBlLO+1){
                    bUpperLO.SetPoint(CounterBlLO,xVal,yVal);
                }
                else{
                    bUpperLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal>yValLowUp){
                        bUpperLO.SetPoint(CounterBlLO,xVal,yVal);
                    }
                }
                if(BestSoFarLO){
                    bBestLO.SetPoint(CounterBlLO,xVal,yVal);
                }
                CounterBlLO++;
            }
            if(Potential==Id_NLO){
                if(bLowerNLO.GetN()<CounterBlNLO+1){
                    bLowerNLO.SetPoint(CounterBlNLO,xVal,yVal);
                }
                else{
                    bLowerNLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal<yValLowUp){
                        bLowerNLO.SetPoint(CounterBlNLO,xVal,yVal);
                    }
                }
                if(bUpperNLO.GetN()<CounterBlNLO+1){
                    bUpperNLO.SetPoint(CounterBlNLO,xVal,yVal);
                }
                else{
                    bUpperNLO.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal>yValLowUp){
                        bUpperNLO.SetPoint(CounterBlNLO,xVal,yVal);
                    }
                }
                if(BestSoFarNLO){
                    bBestNLO.SetPoint(CounterBlNLO,xVal,yVal);
                }
                CounterBlNLO++;
            }
        }
        NumAcceptedEntries++;
        delete FileGraph;
    }

    TH1D* hProjection;
    for(unsigned uPoint=0; uPoint<DevPlotNumBins; uPoint++){
        double kVal = hVarDeviationLO->GetXaxis()->GetBinCenter(uPoint+1);
        hProjection = hVarDeviationLO->ProjectionY("hProjection", uPoint+1,uPoint+1);
        gVarDeviationLO_Stdv.SetPoint(uPoint,kVal,hProjection->GetMean());
        //if( (kVal>265&&kVal<295)||(uPoint==1||uPoint==2||uPoint==3||uPoint==6||uPoint==7||uPoint==8) )
        //{gVarDeviationLO_Stdv.SetPointError(uPoint,0,20);}
        //else if( (uPoint==4||uPoint==5) )
        //{gVarDeviationLO_Stdv.SetPointError(uPoint,0,4);}
        //else
        {gVarDeviationLO_Stdv.SetPointError(uPoint,0,hProjection->GetStdDev());}
        delete hProjection;

        hProjection = hVarDeviationNLO->ProjectionY("hProjection", uPoint+1,uPoint+1);

        gVarDeviationNLO_Stdv.SetPoint(uPoint,kVal,hProjection->GetMean());
        //if( (kVal>265&&kVal<295)||(uPoint==1||uPoint==2||uPoint==3||uPoint==6||uPoint==7||uPoint==8) )
        //{gVarDeviationNLO_Stdv.SetPointError(uPoint,0,20);}
        //else if( (uPoint==4||uPoint==5) )
        //{gVarDeviationNLO_Stdv.SetPointError(uPoint,0,4);}
        //else
        {gVarDeviationNLO_Stdv.SetPointError(uPoint,0,hProjection->GetStdDev());}
        delete hProjection;
    }


    printf("BestNsigma_LO = %.2f (chi2ndf=%.2f)\n",BestNsigma_LO,BestChi2NdfLO);
    printf("BestNsigma_NLO = %.2f (chi2ndf=%.2f)\n",BestNsigma_NLO,BestChi2NdfNLO);

    printf("DefaultChi2NdfLO=%.2f (%.2f)\n",DefaultChi2NdfLO,DefaultNsigmaLO);
    printf("DefaultChi2NdfNLO=%.2f (%.2f)\n",DefaultChi2NdfNLO,DefaultNsigmaNLO);

    printf("MeanChi2NdfLO=%.2f (%.2f)\n",MeanChi2NdfLO,MeanNsigmaLO);
    printf("MeanChi2NdfNLO=%.2f (%.2f)\n",MeanChi2NdfNLO,MeanNsigmaNLO);

    printf("StdvChi2NdfLO=%.2f (%.2f)\n",StdvChi2NdfLO,StdvNsigmaLO);
    printf("StdvChi2NdfNLO=%.2f (%.2f)\n",StdvChi2NdfNLO,StdvNsigmaNLO);

    printf("UniformStdvChi2NdfLO=%.2f (%.2f)\n",UniformStdvChi2NdfLO,UniformStdvNsigmaLO);
    printf("UniformStdvChi2NdfNLO=%.2f (%.2f)\n",UniformStdvChi2NdfNLO,UniformStdvNsigmaNLO);

    printf("MinChi2NdfLO=%.2f (%.2f)\n",MinChi2NdfLO,MinNsigmaLO);
    printf("MinChi2NdfNLO=%.2f (%.2f)\n",MinChi2NdfNLO,MinNsigmaNLO);

    printf("MaxChi2NdfLO=%.2f (%.2f)\n",MaxChi2NdfLO,MaxNsigmaLO);
    printf("MaxChi2NdfNLO=%.2f (%.2f)\n",MaxChi2NdfNLO,MaxNsigmaNLO);

    printf("--- CENTRAL INTERVAL (68%%) ---\n");
    printf(" Chi2NdfLO: %.2f, %.2f, %.2f\n",LowCi_Chi2NdfLO,Median_Chi2NdfLO,UpCi_Chi2NdfLO);
    printf(" Chi2NdfNLO: %.2f, %.2f, %.2f\n",LowCi_Chi2NdfNLO,Median_Chi2NdfNLO,UpCi_Chi2NdfNLO);
    printf(" NsigmaLO: %.2f, %.2f, %.2f\n",LowCi_NsigmaLO,Median_NsigmaLO,UpCi_NsigmaLO);
    printf(" NsigmaNLO: %.2f, %.2f, %.2f\n",LowCi_NsigmaNLO,Median_NsigmaNLO,UpCi_NsigmaNLO);

    TFile* fPlot = new TFile(FitSystFolder+"PLOT/fPlot_"+WhichBaseline+".root","recreate");
    gLowerLO.Write();
    gUpperLO.Write();
    gBestLO.Write();
    gLowerNLO.Write();
    gUpperNLO.Write();
    gBestNLO.Write();
    hChi2NdfLO->Write();
    hChi2NdfNLO->Write();
    hNsigmaLO->Write();
    hNsigmaNLO->Write();
    hVarDeviationLO->Write();
    gVarDeviationLO_Stdv.Write();
    hVarDevRatioLO->Write();
    hVarDeviationNLO->Write();
    gVarDeviationNLO_Stdv.Write();
    hVarDevRatioNLO->Write();

    gStyle->SetCanvasPreferGL(1);
    SetStyle();
    const float right = 0.025;
    const float top = 0.025;

    TGraphErrors *grFemto_LO;
    grFemto_LO = FemtoModelFitBandsSimple(&gLowerLO, &gUpperLO);
    grFemto_LO->SetFillColorAlpha(fColors[3],0.3);
    grFemto_LO->SetLineColor(fColors[3]);
    grFemto_LO->SetLineWidth(5);
    TGraphErrors *grOuterBl_LO;
    grOuterBl_LO = FemtoModelFitBandsSimple(&bOuterLowerLO, &bOuterUpperLO);
    grOuterBl_LO->SetFillColorAlpha(fColors[3],0.3);
    grOuterBl_LO->SetLineColor(fColors[3]);
    grOuterBl_LO->SetLineWidth(5);
    gBestLO.SetLineColor(fColors[3]);
    gBestLO.SetLineWidth(5);
    bOuterBestLO.SetLineColor(fColors[3]);
    bOuterBestLO.SetLineWidth(5);
    bBestLO.SetLineColor(fColors[3]);
    bBestLO.SetLineWidth(4);
    bBestLO.SetLineStyle(2);

    TGraphErrors *grFemto_NLO;
    grFemto_NLO = FemtoModelFitBandsSimple(&gLowerNLO, &gUpperNLO);
    grFemto_NLO->SetFillColorAlpha(fColors[1],0.3);
    grFemto_NLO->SetLineColor(fColors[1]);
    grFemto_NLO->SetLineWidth(5);
    TGraphErrors *grOuterBl_NLO;
    grOuterBl_NLO = FemtoModelFitBandsSimple(&bOuterLowerNLO, &bOuterUpperNLO);
    grOuterBl_NLO->SetFillColorAlpha(fColors[1],0.3);
    grOuterBl_NLO->SetLineColor(fColors[1]);
    grOuterBl_NLO->SetLineWidth(5);
    gBestNLO.SetLineColor(fColors[1]);
    gBestNLO.SetLineWidth(5);
    bOuterBestNLO.SetLineColor(fColors[1]);
    bOuterBestNLO.SetLineWidth(5);
    bBestNLO.SetLineColor(fColors[1]);
    bBestNLO.SetLineWidth(4);
    bBestNLO.SetLineStyle(7);

    DLM_SubPads DlmPad(720,720);
    DlmPad.AddSubPad(0,1,0.33,1);
    DlmPad.AddSubPad(0,1,0,0.33);
    DlmPad.SetMargin(0,0.12,0.02,0.0,0.02);
    DlmPad.SetMargin(1,0.12,0.02,0.09,0.0);
    DlmPad.cd(0);

    //TCanvas *Can_CF_pL = new TCanvas("pL","pL", 0,0,720,720);
    //Can_CF_pL->SetRightMargin(right);
    //Can_CF_pL->SetTopMargin(top);

    hData->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hData->GetXaxis()->SetRangeUser(0, 320);
    hData->GetXaxis()->SetNdivisions(505);
    hData->GetYaxis()->SetRangeUser(0.85, 2.3);
    hData->SetFillColor(fFillColors[0]);
    SetStyleHisto2(hData,2,0);
    //hData->GetYaxis()->SetTitleOffset(1.0);
    hData->Draw();

    TGraphErrors *Tgraph_syserror = DrawSystematicError_FAST(hData, hSyst, fSyst, 3);
    Tgraph_syserror->SetLineColor(kWhite);

    //baselineLL->Draw("same");

    if(grFemto_LO) {grFemto_LO->Draw("3 same");}
    if(grOuterBl_LO&&MODE==1) {grOuterBl_LO->Draw("3 same");}
    if(grFemto_NLO) grFemto_NLO->Draw("3 same");
    if(grOuterBl_NLO&&MODE==1) {grOuterBl_NLO->Draw("3 same");}
    if(grFemto_LO) gBestLO.Draw("l same");
    if(grFemto_LO&&MODE==1) {bOuterBestLO.Draw("l same");}
    if(grFemto_NLO) gBestNLO.Draw("l same");
    if(grFemto_NLO&&MODE==1) {bOuterBestNLO.Draw("l same");}
    if(grFemto_LO&&MODE==0) {bBestLO.Draw("l same");}
    if(grFemto_NLO&&MODE==0) {bBestNLO.Draw("l same");}
    hData->Draw("same");

    Tgraph_syserror->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror->Draw("2 same");
    //hData->Draw("pe same");

    unsigned NumRows=4;
    TLegend *legend = new TLegend(0.49,0.72-0.05*NumRows,0.73,0.72);//lbrt
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(gStyle->GetTextSize()*0.90);
    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hData->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    legend->AddEntry(hCk_Fake, "p#minus#Lambda #oplus #bar{p}#minus#bar{#Lambda} pairs", "fpe");
//  legend->AddEntry(hist_CF_LL_ALAL_exp[2], "with Syst. uncertainties", "");
    //legend->AddEntry(baselineLL,"Baseline","l");
    legend->AddEntry(grFemto_LO,"Femtoscopic fit (#chiEFT LO)","l");
    legend->AddEntry(grFemto_NLO,"Femtoscopic fit (#chiEFT NLO)","l");
    legend->Draw("same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
    //if(WhichDataSet==0) BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 13 TeV"));
    //else if(WhichDataSet==1) BeamText.DrawLatex(0.55, 0.825, Form("p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    //else BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 7 TeV"));
    BeamText.DrawLatex(0.50, 0.915, "ALICE Preliminary");
    //if(DataSample=="pp13TeV_MB_Run2paper") BeamText.DrawLatex(0.50, 0.86, "pp #sqrt{#it{s}} = 13 TeV");
    //else if(DataSample.Contains("pPb")&&DataSample.Contains("5TeV")) BeamText.DrawLatex(0.50, 0.86, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    //else if(DataSample=="pp13TeV_HM_March19") BeamText.DrawLatex(0.50, 0.86, "pp (HM) #sqrt{#it{s}} = 13 TeV");
    //else BeamText.DrawLatex(0.50, 0.86, "ALICE pp #sqrt{#it{s}} = 7 TeV");
    BeamText.DrawLatex(0.50, 0.860, "pp (HM) #sqrt{#it{s}} = 13 TeV");

    TLatex BeamTextSource;
    BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamTextSource.SetNDC(kTRUE);
    BeamTextSource.DrawLatex(0.50, 0.805, LegendSource_line1);
    BeamTextSource.DrawLatex(0.50, 0.750, LegendSource_line2);

//INLET -------------------------------------------------------------------------------------------------------------------

    TH1F* DataHisto_Inlet = (TH1F*)hData->Clone("DataHisto_Inlet");
    DataHisto_Inlet->SetMarkerSize(hData->GetMarkerSize()*0.67);
    DataHisto_Inlet->SetLineWidth(hData->GetLineWidth()*0.67);
    DataHisto_Inlet->GetXaxis()->SetTitleSize(hData->GetXaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetLabelSize(hData->GetXaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetRangeUser(120, MODE==1?600:560);
    DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

    DataHisto_Inlet->GetYaxis()->SetTitleSize(hData->GetYaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetLabelSize(hData->GetYaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetTitleOffset(hData->GetYaxis()->GetTitleOffset()*0.67);
    DataHisto_Inlet->GetYaxis()->SetRangeUser(0.962, 1.048);

    TGraph* grFemto_LO_Inlet = (TGraph*)grFemto_LO->Clone("grFemto_LO_Inlet");
    grFemto_LO_Inlet->SetLineWidth(grFemto_LO->GetLineWidth()*0.67);

    TGraph* grFemto_NLO_Inlet = (TGraph*)grFemto_NLO->Clone("grFemto_NLO_Inlet");
    grFemto_NLO_Inlet->SetLineWidth(grFemto_NLO->GetLineWidth()*0.67);

    TGraph* gBestLO_Inlet = (TGraph*)gBestLO.Clone("gBestLO_Inlet");
    gBestLO_Inlet->SetLineWidth(gBestLO.GetLineWidth()*0.67);

    TGraph* gBestNLO_Inlet = (TGraph*)gBestNLO.Clone("gBestNLO_Inlet");
    gBestNLO_Inlet->SetLineWidth(gBestNLO.GetLineWidth()*0.67);

    const double fXMinInlet=0.30;
    const double fYMinInlet=0.12;
    const double fXMaxInlet=0.95;
    const double fYMaxInlet=0.50;
    TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                             fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    DataHisto_Inlet->Draw();
    //if(grFemto_LO) grFemto_LO_Inlet->Draw("3 same");
    //if(grFemto_NLO) grFemto_NLO_Inlet->Draw("3 same");
    //if(grFemto_LO) gBestLO_Inlet->Draw("l same");
    //if(grFemto_NLO) gBestNLO_Inlet->Draw("l same");
    if(grFemto_LO) {grFemto_LO->Draw("3 same");}
    if(grOuterBl_LO&&MODE==1) {grOuterBl_LO->Draw("3 same");}
    if(grFemto_NLO) grFemto_NLO->Draw("3 same");
    if(grOuterBl_NLO&&MODE==1) {grOuterBl_NLO->Draw("3 same");}
    if(grFemto_LO) gBestLO.Draw("l same");
    if(grFemto_LO&&MODE==1) {bOuterBestLO.Draw("l same");}
    if(grFemto_NLO) gBestNLO.Draw("l same");
    if(grFemto_NLO&&MODE==1) {bOuterBestNLO.Draw("l same");}
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror->Draw("2 same");
    if(grFemto_LO&&MODE==0) {bBestLO.Draw("l same");}
    if(grFemto_NLO&&MODE==0) {bBestNLO.Draw("l same");}

    DlmPad.cd(1);
    TH1F* hAxis = new TH1F("hAxis", "hAxis", 12, 0, 324);
    hAxis->SetStats(false);
    hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");

    hAxis->GetYaxis()->SetRangeUser(-11.5, 11.5);
    //hData->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    //hData->GetXaxis()->SetRangeUser(0, 320);
    hAxis->GetXaxis()->SetNdivisions(505);
    //hData->GetYaxis()->SetRangeUser(0.85, 2.3);
    //hData->SetFillColor(fFillColors[0]);
    SetStyleHisto2(hAxis,2,0,2);
    //hData->GetYaxis()->SetTitleOffset(1.0);
    hAxis->Draw("");

    gVarDeviationLO_Stdv.SetFillColorAlpha(kGreen+3,0.75);
    gVarDeviationLO_Stdv.Draw("3,same");
    gVarDeviationNLO_Stdv.SetFillColorAlpha(kRed+1,0.75);
    gVarDeviationNLO_Stdv.Draw("3,same");

    DlmPad.GetCanvas()->SaveAs(FitSystFolder+TString::Format("PLOT/Can_CF_pL_%i.pdf",WhichBaseline));

    delete hVarDeviationLO;
    delete hVarDeviationNLO;
    delete hVarDevRatioLO;
    delete hVarDevRatioNLO;
    delete legend;
    delete hCk_Fake;
    delete Tgraph_syserror;
    delete grFemto_LO;
    delete grFemto_NLO;
    delete grOuterBl_LO;
    delete grOuterBl_NLO;
    //delete Can_CF_pL;
    delete hChi2NdfLO;
    delete hChi2NdfNLO;
    delete hNsigmaLO;
    delete hNsigmaNLO;
    delete fPlot;
    delete DataFile;
    delete SystFile;
    delete ntFile;
/*
    //const int NumConfigs = 2;
    //int WhichConfiguration[NumConfigs] = {0,2};

    //const int NumConfigs = 3;
    //int WhichConfiguration[NumConfigs] = {0,1,2};

    TString FitResultHistoName1 = "pL_AV18_"+pL_Model1+"_pXim_HALQCD1";
    TString FitResultHistoName2 = "pL_AV18_"+pL_Model2+"_pXim_HALQCD1";

    TFile* DataFile = new TFile(DataFileName,"read");
    if(!DataFile){
        printf("Problem with %s\n",DataFileName.Data());
        return;
    }
    TH1F* DataHisto = (TH1F*)DataFile->Get(DataHistoName);
    if(!DataHisto){
        printf("Problem with %s\n",DataHistoName.Data());
        return;
    }
    TGraph* FitResult1 = (TGraph*)DataFile->Get(FitResultHistoName1);
    TGraph* FitResult2 = (TGraph*)DataFile->Get(FitResultHistoName2);

    TGraph gDownLimit;
    TGraph gUpLimit;

    TFile* DataSystFile = new TFile(DataSystFileName);
    TH1F* DataSystHisto = (TH1F*)DataSystFile->Get(DataSystHistoName);
    TF1* ExtrapolOLD = NULL;
    if(!DataSystHisto){
        ExtrapolOLD = (TF1*)DataSystFile->Get("RelSysPLUnbinned");
    }
//printf("DataSystFile=%p\n",DataSystFile);
//printf("DataSystHisto=%p\n",DataSystHisto);
//printf("ExtrapolOLD=%p\n",ExtrapolOLD);

    gStyle->SetCanvasPreferGL(1);
    SetStyle();
    const float right = 0.025;
    const float top = 0.025;

    TGraphErrors *grFemtoLL_1;
    grFemtoLL_1 = FemtoModelFitBandsSimple(FitResult1, FitResult1);
    grFemtoLL_1->SetFillColor(fColors[3]);
    grFemtoLL_1->SetLineColor(fColors[3]);
    grFemtoLL_1->SetLineWidth(5);

    TGraphErrors *grFemtoLL_2;
    grFemtoLL_2 = FemtoModelFitBandsSimple(FitResult2, FitResult2);
    grFemtoLL_2->SetFillColor(fColors[1]);
    grFemtoLL_2->SetLineColor(fColors[1]);
    grFemtoLL_2->SetLineWidth(5);

    TCanvas *Can_CF_LL = new TCanvas("LL","LL", 0,0,650,550);
    Can_CF_LL->SetRightMargin(right);
    Can_CF_LL->SetTopMargin(top);

    TF1 *baselineLL = new TF1("baselineLL", "pol1", 0, 1000);
    baselineLL->SetParameter(0,parA);
    baselineLL->SetParameter(1,parB);
    baselineLL->SetLineStyle(2);
    baselineLL->SetLineColor(fFillColors[6]);
    baselineLL->SetLineWidth(4);

    TGraph *grFakeLL = new TGraph();
    grFakeLL->SetLineColor(fColors[5]);
    grFakeLL->SetLineWidth(4);

    unsigned NumMomBins = DataHisto->GetNbinsX();
    double MinMom = DataHisto->GetBinLowEdge(1);
    double MaxMom = DataHisto->GetXaxis()->GetBinUpEdge(NumMomBins);
    TH1F* hGeV = new TH1F("hGeV","hGeV",NumMomBins,MinMom*0.001,MaxMom*0.001);
    for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
        hGeV->SetBinContent(uBin,DataHisto->GetBinContent(uBin));
        hGeV->SetBinError(uBin,DataHisto->GetBinError(uBin));
    }
    hGeV->SetTitle("; #it{k*} (GeV/#it{c}); #it{C}(#it{k*})");
    hGeV->GetXaxis()->SetRangeUser(0, 0.32);
    hGeV->GetXaxis()->SetNdivisions(505);
    hGeV->GetYaxis()->SetRangeUser(0.85, 2.3);
    hGeV->SetFillColor(fFillColors[0]);
    SetStyleHisto(hGeV,2,0);
    hGeV->SetMarkerSize(1.5);
    //hGeV->Draw();

    DataHisto->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    DataHisto->GetXaxis()->SetRangeUser(0, 320);
    DataHisto->GetXaxis()->SetNdivisions(505);
    DataHisto->GetYaxis()->SetRangeUser(0.85, 2.3);
    DataHisto->SetFillColor(fFillColors[0]);
    SetStyleHisto(DataHisto,2,0);
    DataHisto->Draw();

    TGraphErrors *Tgraph_syserror_LL_ALAL = DrawSystematicError_FAST(DataHisto, DataSystHisto, ExtrapolOLD, 3);
    Tgraph_syserror_LL_ALAL->SetLineColor(kWhite);

    //baselineLL->Draw("same");

    if(grFemtoLL_1) grFemtoLL_1->Draw("l3 same");
    if(grFemtoLL_2) grFemtoLL_2->Draw("l3 same");
    //hGeV->Draw("same");
    DataHisto->Draw("same");

    Tgraph_syserror_LL_ALAL->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror_LL_ALAL->Draw("2 same");
    //DataHisto->Draw("pe same");

    unsigned NumRows=4;
    TLegend *legLL2 = new TLegend(0.49,0.76-0.04*NumRows,0.73,0.76);//lbrt
    legLL2->SetBorderSize(0);
    legLL2->SetTextFont(42);
    legLL2->SetTextSize(gStyle->GetTextSize()*0.70);
    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hGeV->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    legLL2->AddEntry(hCk_Fake, "p#minus#Lambda #oplus #bar{p}#minus#bar{#Lambda} pairs", "fpe");
//  legLL2->AddEntry(hist_CF_LL_ALAL_exp[2], "with Syst. uncertainties", "");
    //legLL2->AddEntry(baselineLL,"Baseline","l");
    legLL2->AddEntry(grFemtoLL_1,"Femtoscopic fit (#chiEFT LO)","l");
    legLL2->AddEntry(grFemtoLL_2,"Femtoscopic fit (#chiEFT NLO)","l");
    legLL2->Draw("same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.80);
    BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
    //if(WhichDataSet==0) BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 13 TeV"));
    //else if(WhichDataSet==1) BeamText.DrawLatex(0.55, 0.825, Form("p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    //else BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 7 TeV"));
    BeamText.DrawLatex(0.50, 0.91, "ALICE Preliminary");
    if(DataSample=="pp13TeV_MB_Run2paper") BeamText.DrawLatex(0.50, 0.86, "pp #sqrt{#it{s}} = 13 TeV");
    else if(DataSample.Contains("pPb")&&DataSample.Contains("5TeV")) BeamText.DrawLatex(0.50, 0.86, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    else if(DataSample=="pp13TeV_HM_March19") BeamText.DrawLatex(0.50, 0.86, "pp (HM) #sqrt{#it{s}} = 13 TeV");
    else BeamText.DrawLatex(0.50, 0.86, "ALICE pp #sqrt{#it{s}} = 7 TeV");

    TLatex BeamTextSource;
    BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.70);
    BeamTextSource.SetNDC(kTRUE);
    BeamTextSource.DrawLatex(0.50, 0.81, LegendSource_line1);



//INLET -------------------------------------------------------------------------------------------------------------------

    TH1F* DataHisto_Inlet = (TH1F*)DataHisto->Clone("DataHisto_Inlet");
    DataHisto_Inlet->SetMarkerSize(DataHisto->GetMarkerSize()*0.67);
    DataHisto_Inlet->SetLineWidth(DataHisto->GetLineWidth()*0.67);
    DataHisto_Inlet->GetXaxis()->SetTitleSize(DataHisto->GetXaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetLabelSize(DataHisto->GetXaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetRangeUser(120, 320);
    DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

    DataHisto_Inlet->GetYaxis()->SetTitleSize(DataHisto->GetYaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetLabelSize(DataHisto->GetYaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetTitleOffset(DataHisto->GetYaxis()->GetTitleOffset()*0.67);
    DataHisto_Inlet->GetYaxis()->SetRangeUser(0.98, 1.045);

    TGraph* grFemtoLL_1_Inlet = (TGraph*)grFemtoLL_1->Clone("grFemtoLL_1_Inlet");
    grFemtoLL_1_Inlet->SetLineWidth(grFemtoLL_1->GetLineWidth()*0.67);

    TGraph* grFemtoLL_2_Inlet = (TGraph*)grFemtoLL_2->Clone("grFemtoLL_2_Inlet");
    grFemtoLL_2_Inlet->SetLineWidth(grFemtoLL_2->GetLineWidth()*0.67);

    const double fXMinInlet=0.35;
    const double fYMinInlet=0.25;
    const double fXMaxInlet=0.95;
    const double fYMaxInlet=0.57;
    TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                             fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    DataHisto_Inlet->Draw();
    if(grFemtoLL_1) grFemtoLL_1_Inlet->Draw("l3 same");
    if(grFemtoLL_2) grFemtoLL_2_Inlet->Draw("l3 same");
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror_LL_ALAL->Draw("2 same");

    Can_CF_LL->Print(WorkFolder+TString::Format("PlotLamLamFit_%s_%s_%s.pdf",DataSample.Data(),SourceType.Data(),FittingMode_pp.Data()));

    delete inset_pad;
    delete DataHisto_Inlet;
    delete grFemtoLL_1_Inlet;
    delete grFemtoLL_2_Inlet;

    delete grFemtoLL_1;
    delete grFemtoLL_2;
    delete Can_CF_LL;
    delete baselineLL;
    delete legLL2;
    delete grFakeLL;
    delete hGeV;
    delete hCk_Fake;

    delete DataFile;
    delete DataSystFile;
*/
}

void Plot_pp_FASTsyst( const int& WhichBaseline ){

//MODE = 0 -> the coupled fit
//MODE = 1 -> the one without coupled
    int MODE = 0;
    double MinMomentum=0;
    double MaxMomentum=336;
    unsigned DefaultIterAV18=0;

    //we take 10HM not to redo further plots
    TString DataFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/CFOutput_pp.root";
    //TString DataFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample12HM/CFOutput_pp_0.root";
    TString DataHistoName = "hCk_ReweightedMeV_0";
    TFile* DataFile = new TFile(DataFileName,"read");
    TH1F* hData = (TH1F*)DataFile->Get(DataHistoName);

    TString DataSystFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample12HM/Systematics_pp.root";
    TString DataSystHistoName = "hRelSyst";//"RelSysPLUnbinned" "hRelSyst" "SystErrRel"
    TFile* SystFile = new TFile(DataSystFileName,"read");
    TH1F* hSyst = (TH1F*)SystFile->Get(DataSystHistoName);
    TF1* fSyst = (TF1*)SystFile->Get(DataSystHistoName);
//SystematicsAdd_100419_2 pm15
//SystematicsAdd_120419 pm10
//SystematicsAddMeson_230419
//SystematicsAdd_310519
//SystematicsAdd_010519
    TString FitSystFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp/SystematicsDef_130519/";
    TString FitSystFileName = FitSystFolder+"NTfile_0.root";

    TFile* ntFile = new TFile(FitSystFileName,"read");
    TNtuple* ntResult = (TNtuple*)ntFile->Get("ntResult");
    unsigned NumNtEntries = ntResult->GetEntries();

    Float_t Iter;
    Float_t Config;
    Float_t Data;
    Float_t SourceType;
    Float_t SourceScale;
    Float_t Potential;
    Float_t Baseline;
    Float_t FemRan;
    Float_t FitRan;
    Float_t pFrac;
    Float_t LamFrac;
    Float_t kcVar;
    Float_t mTbin;
    Float_t FemtoMin;
    Float_t FemtoMax;
    Float_t BlMin;
    Float_t BlMax;
    Float_t p_a;
    Float_t e_a;
    Float_t p_b;
    Float_t e_b;
    Float_t p_c;
    Float_t e_c;
    Float_t p_Cl;
    Float_t e_Cl;
    Float_t p_kc;
    Float_t e_kc;
    Float_t p_sor0;
    Float_t e_sor0;
    Float_t p_sor1;
    Float_t e_sor1;
    Float_t chi2;
    Float_t ndf;
    Float_t chi2_352;
    Float_t ndf_352;

    ntResult->SetBranchAddress("Iter",&Iter);
    ntResult->SetBranchAddress("Config",&Config);
    ntResult->SetBranchAddress("Data",&Data);
    ntResult->SetBranchAddress("SourceType",&SourceType);
    ntResult->SetBranchAddress("SourceScale",&SourceScale);
    ntResult->SetBranchAddress("Potential",&Potential);
    ntResult->SetBranchAddress("Baseline",&Baseline);
    ntResult->SetBranchAddress("FemRan",&FemRan);
    ntResult->SetBranchAddress("FitRan",&FitRan);
    ntResult->SetBranchAddress("pFrac",&pFrac);
    ntResult->SetBranchAddress("LamFrac",&LamFrac);
    ntResult->SetBranchAddress("kcVar",&kcVar);
    ntResult->SetBranchAddress("mTbin",&mTbin);
    ntResult->SetBranchAddress("FemtoMin",&FemtoMin);
    ntResult->SetBranchAddress("FemtoMax",&FemtoMax);
    ntResult->SetBranchAddress("BlMin",&BlMin);
    ntResult->SetBranchAddress("BlMax",&BlMax);
    ntResult->SetBranchAddress("p_a",&p_a);
    ntResult->SetBranchAddress("e_a",&e_a);
    ntResult->SetBranchAddress("p_b",&p_b);
    ntResult->SetBranchAddress("e_b",&e_b);
    ntResult->SetBranchAddress("p_c",&p_c);
    ntResult->SetBranchAddress("e_c",&e_c);
    ntResult->SetBranchAddress("p_Cl",&p_Cl);
    ntResult->SetBranchAddress("e_Cl",&e_Cl);
    ntResult->SetBranchAddress("p_kc",&p_kc);
    ntResult->SetBranchAddress("e_kc",&e_kc);
    ntResult->SetBranchAddress("p_sor0",&p_sor0);
    ntResult->SetBranchAddress("e_sor0",&e_sor0);
    ntResult->SetBranchAddress("p_sor1",&p_sor1);
    ntResult->SetBranchAddress("e_sor1",&e_sor1);
    ntResult->SetBranchAddress("chi2",&chi2);
    ntResult->SetBranchAddress("ndf",&ndf);
    ntResult->SetBranchAddress("chi2_352",&chi2_352);
    ntResult->SetBranchAddress("ndf_352",&ndf_352);

    TGraph gLowerAV18;
    gLowerAV18.SetName("gLowerAV18");
    TGraph gUpperAV18;
    gUpperAV18.SetName("gUpperAV18");
    TGraph gBestAV18;
    gBestAV18.SetName("gBestAV18");

    //baseline
    TGraph bOuterLowerAV18;
    bOuterLowerAV18.SetName("bOuterLowerAV18");
    TGraph bOuterUpperAV18;
    bOuterUpperAV18.SetName("bOuterUpperAV18");
    TGraph bOuterBestAV18;
    bOuterBestAV18.SetName("bOuterBestAV18");
    TGraph bLowerAV18;
    bLowerAV18.SetName("bLowerAV18");
    TGraph bUpperAV18;
    bUpperAV18.SetName("bUpperAV18");
    TGraph bBestAV18;
    bBestAV18.SetName("bBestAV18");

    double BestChi2NdfAV18 = 1e6;
    double WorstChi2NdfAV18 = 0;
    unsigned NumAcceptedEntries = 0;

    double GlobalBestPval_AV18 = 0;



    TH1F* hChi2NdfAV18 = new TH1F("hChi2NdfAV18","hChi2NdfAV18",1024,0,64);
    TH1F* hNsigmaAV18 = new TH1F("hNsigmaAV18","hNsigmaAV18",1024,0,64);

    //bin-by-bin deviation in nSigma
    //we use TH2F to save all variations
    //than we take the mean and stdv to define upper and lower limit
    const unsigned DevPlotNumBins = 150;
    TH2F* hVarDeviationAV18 = new TH2F("hVarDeviationAV18","hVarDeviationAV18",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,1024,-40,40);
    TH2F* hVarDevRatioAV18 = new TH2F("hVarDevRatioAV18","hVarDevRatioAV18",DevPlotNumBins,0,hData->GetBinWidth(1)*DevPlotNumBins,128,0.85,1.15);

    TGraphErrors gVarDeviationAV18_Stdv;
    gVarDeviationAV18_Stdv.SetName("gVarDeviationAV18_Stdv");
    gVarDeviationAV18_Stdv.SetLineWidth(0);
    gVarDeviationAV18_Stdv.SetLineColor(kBlue+2);
    gVarDeviationAV18_Stdv.SetMarkerSize(0);
    gVarDeviationAV18_Stdv.SetMarkerStyle(20);
    gVarDeviationAV18_Stdv.SetMarkerColor(kBlue+2);
    gVarDeviationAV18_Stdv.SetFillColor(kBlue+2);


    double DefaultChi2NdfAV18=0;
    double DefaultNsigmaAV18=0;

    for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
        ntFile->cd();
        ntResult->GetEntry(uEntry);
        int ITER = int(Iter);
        int CONFIG = int(Config);
        int Id_AV18;
        if(CONFIG<10){
            Id_AV18 = 0;
        }
        else if(CONFIG<20){
            Id_AV18 = 0;
        }

//printf("Iter=%i; CONFIG=%i; LO=%i; NLO=%i\n",ITER,CONFIG,Id_LO,Id_NLO);
        //! Conditions

        //if(CONFIG!=0&&CONFIG!=2) continue; //select Gauss+Reso and Gauss
        if(CONFIG!=0) continue; //select Gauss+Reso
        //if(CONFIG!=1) continue; //select Levy+Reso
        //if(CONFIG!=2) continue; //select Gauss
        if(round(Baseline)!=WhichBaseline) continue;
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        if(Id_AV18==Potential){
            hChi2NdfAV18->Fill(chi2_352/ndf_352);
            double nSigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(chi2_352,int(ndf_352)));
            if(nSigma==0) nSigma = 100;
            hNsigmaAV18->Fill(nSigma);
            if(DefaultIterAV18==int(Iter)){
                DefaultChi2NdfAV18 = chi2_352/ndf_352;
                DefaultNsigmaAV18 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(chi2_352,int(ndf_352)));
            }
        }
        else{
            printf("PROBLEM\n");
        }
    }

    double MeanChi2NdfAV18 = hChi2NdfAV18->GetMean();
    double StdvChi2NdfAV18 = hChi2NdfAV18->GetStdDev();
    double LenChi2NdfAV18;
    double UniformStdvChi2NdfAV18;
    double MinChi2NdfAV18=1e6;
    double MaxChi2NdfAV18=-1e6;

    double MeanNsigmaAV18 = hNsigmaAV18->GetMean();
    double StdvNsigmaAV18 = hNsigmaAV18->GetStdDev();
    double LenNsigmaAV18;
    double UniformStdvNsigmaAV18;
    double MinNsigmaAV18=1e6;
    double MaxNsigmaAV18=-1e6;


    for(unsigned uBin=0; uBin<hChi2NdfAV18->GetNbinsX(); uBin++){
        if(hChi2NdfAV18->GetBinContent(uBin+1)&&MinChi2NdfAV18>hChi2NdfAV18->GetBinCenter(uBin+1)){MinChi2NdfAV18=hChi2NdfAV18->GetBinCenter(uBin+1);}
        if(hChi2NdfAV18->GetBinContent(uBin+1)&&MaxChi2NdfAV18<hChi2NdfAV18->GetBinCenter(uBin+1)){MaxChi2NdfAV18=hChi2NdfAV18->GetBinCenter(uBin+1);}
    }
    LenChi2NdfAV18 = MaxChi2NdfAV18-MinChi2NdfAV18;
    UniformStdvChi2NdfAV18 = LenChi2NdfAV18/sqrt(12.);

    for(unsigned uBin=0; uBin<hNsigmaAV18->GetNbinsX(); uBin++){
        if(hNsigmaAV18->GetBinContent(uBin+1)&&MinNsigmaAV18>hNsigmaAV18->GetBinCenter(uBin+1)){MinNsigmaAV18=hNsigmaAV18->GetBinCenter(uBin+1);}
        if(hNsigmaAV18->GetBinContent(uBin+1)&&MaxNsigmaAV18<hNsigmaAV18->GetBinCenter(uBin+1)){MaxNsigmaAV18=hNsigmaAV18->GetBinCenter(uBin+1);}
    }
    LenNsigmaAV18 = MaxNsigmaAV18-MinNsigmaAV18;
    UniformStdvNsigmaAV18 = LenNsigmaAV18/sqrt(12.);


    double BestNsigma_AV18 = sqrt(2)*TMath::ErfcInverse(GlobalBestPval_AV18);

    double LowCi_Chi2NdfAV18;
    double UpCi_Chi2NdfAV18;
    double Median_Chi2NdfAV18 = GetCentralInterval(*hChi2NdfAV18,0.68,LowCi_Chi2NdfAV18,UpCi_Chi2NdfAV18,true);

    double LowCi_NsigmaAV18;
    double UpCi_NsigmaAV18;
    double Median_NsigmaAV18 = GetCentralInterval(*hNsigmaAV18,0.68,LowCi_NsigmaAV18,UpCi_NsigmaAV18,true);

    //THANK YOU ROOT, THANKS!!!
    //so if we read the same TNtuple twice, the second iteration is SUPER slow (like c.a. x1000 times).
    //to fix this: delete the TNtuple and open it again
    delete ntResult;
    ntResult = (TNtuple*)ntFile->Get("ntResult");

    ntResult->SetBranchAddress("Iter",&Iter);
    ntResult->SetBranchAddress("Config",&Config);
    ntResult->SetBranchAddress("Data",&Data);
    ntResult->SetBranchAddress("SourceType",&SourceType);
    ntResult->SetBranchAddress("SourceScale",&SourceScale);
    ntResult->SetBranchAddress("Potential",&Potential);
    ntResult->SetBranchAddress("Baseline",&Baseline);
    ntResult->SetBranchAddress("FemRan",&FemRan);
    ntResult->SetBranchAddress("FitRan",&FitRan);
    ntResult->SetBranchAddress("pFrac",&pFrac);
    ntResult->SetBranchAddress("LamFrac",&LamFrac);
    ntResult->SetBranchAddress("kcVar",&kcVar);
    ntResult->SetBranchAddress("mTbin",&mTbin);
    ntResult->SetBranchAddress("FemtoMin",&FemtoMin);
    ntResult->SetBranchAddress("FemtoMax",&FemtoMax);
    ntResult->SetBranchAddress("BlMin",&BlMin);
    ntResult->SetBranchAddress("BlMax",&BlMax);
    ntResult->SetBranchAddress("p_a",&p_a);
    ntResult->SetBranchAddress("e_a",&e_a);
    ntResult->SetBranchAddress("p_b",&p_b);
    ntResult->SetBranchAddress("e_b",&e_b);
    ntResult->SetBranchAddress("p_c",&p_c);
    ntResult->SetBranchAddress("e_c",&e_c);
    ntResult->SetBranchAddress("p_Cl",&p_Cl);
    ntResult->SetBranchAddress("e_Cl",&e_Cl);
    ntResult->SetBranchAddress("p_kc",&p_kc);
    ntResult->SetBranchAddress("e_kc",&e_kc);
    ntResult->SetBranchAddress("p_sor0",&p_sor0);
    ntResult->SetBranchAddress("e_sor0",&e_sor0);
    ntResult->SetBranchAddress("p_sor1",&p_sor1);
    ntResult->SetBranchAddress("e_sor1",&e_sor1);
    ntResult->SetBranchAddress("chi2",&chi2);
    ntResult->SetBranchAddress("ndf",&ndf);
    ntResult->SetBranchAddress("chi2_352",&chi2_352);//chi2<=312 MeV
    ntResult->SetBranchAddress("ndf_352",&ndf_352);//number of data points up to 312 MeV

//TRandom3 rangen(11);
//printf("NumNtEntries=%u\n",NumNtEntries);
    for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
        //if(uEntry%100==0) printf("uEntry=%u\n",uEntry);
//if(rangen.Uniform()>0.01) continue;
        ntFile->cd();
        ntResult->GetEntry(uEntry);

        int ITER = int(Iter);
        int CONFIG = int(Config);
        int Id_AV18;
        if(CONFIG<10){
            Id_AV18 = 0;
        }
        else if(CONFIG<20){
            Id_AV18 = 0;
        }
//printf("Iter=%i; CONFIG=%i; LO=%i; NLO=%i\n",ITER,CONFIG,Id_LO,Id_NLO);
        //! Conditions
//if(rangen.Uniform()>0.01&&ITER!=DefaultIterLO&&ITER!=DefaultIterNLO) continue;
        //if(CONFIG!=0&&CONFIG!=2) continue; //select Gauss+Reso and Gauss
        if(CONFIG!=0) continue; //select Gauss+Reso
        //if(CONFIG!=1) continue; //select Levy+Reso
        //if(CONFIG!=2) continue; //select Gauss
        if(round(Baseline)!=WhichBaseline) continue;
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        if(Potential==Id_AV18&&chi2_352/double(ndf_352)<LowCi_Chi2NdfAV18&&chi2_352/double(ndf_352)>UpCi_Chi2NdfAV18) continue;
        //if(Data!=0) continue;

        bool BestSoFarAV18=false;

        //we want to plot something that makes sense (i.e. is fixed to the data we plot)
        if(Data==0){
            if(Potential==Id_AV18&&(chi2_352/ndf_352)>WorstChi2NdfAV18&&pFrac==0&&LamFrac==0){
                WorstChi2NdfAV18=chi2_352/ndf_352;
            }
            if(Potential==Id_AV18&&(chi2_352/ndf_352)<BestChi2NdfAV18&&pFrac==0&&LamFrac==0){
                BestChi2NdfAV18=chi2_352/ndf_352;
                BestSoFarAV18=true;
            }
            if(Potential==Id_AV18&&TMath::Prob(chi2_352,ndf_352)>GlobalBestPval_AV18){
                GlobalBestPval_AV18=TMath::Prob(chi2_352,ndf_352);
    //printf("GlobalBestPval_AV18=%f\n",GlobalBestPval_AV18);
                //BestSoFarAV18=true;
            }

        }

        TFile* FileGraph = new TFile(TString::Format("%sConfig%i_Iter%i.root",FitSystFolder.Data(),CONFIG,ITER),"read");
        TH1F* hDataVar = (TH1F*)FileGraph->Get("hCk_ReweightedMeV_0");
        TGraph* FitResult = (TGraph*)FileGraph->Get("FitResult_pp");
        TF1* FitBaseline = (TF1*)FileGraph->Get("fBaseline");

        double CkMaxFit;
        double kMaxFit;
        FitResult->GetPoint(FitResult->GetN()-1,kMaxFit,CkMaxFit);

        for(unsigned uBin=0; uBin<DevPlotNumBins; uBin++){
            double kValue=hVarDeviationAV18->GetXaxis()->GetBinCenter(uBin+1);
            double dataValue = hDataVar->GetBinContent(uBin+1);
            double fitValue = FitResult->Eval(kValue);
            if(kValue>kMaxFit){
                fitValue = CkMaxFit;
            }
            double dataError = hDataVar->GetBinError(uBin+1);
            double Nsigma=(fitValue-dataValue)/(dataError);
            double Ratio=dataValue/fitValue;

            if(Potential==Id_AV18){
                hVarDeviationAV18->Fill(kValue,Nsigma);
                hVarDevRatioAV18->Fill(kValue,Ratio);
            }
        }

        double xVal,xValLowUp;
        double yVal,yValLowUp;
        TGraph& gLower = gLowerAV18;
        TGraph& gUpper = gUpperAV18;

        //TGraph& gDefault = Potential==1?gDefaultLO:gDefaultNLO;
        for(int iPoint=0; iPoint<FitResult->GetN(); iPoint++){
            FitResult->GetPoint(iPoint,xVal,yVal);
            if(gLower.GetN()<iPoint+1){
                gLower.SetPoint(iPoint,xVal,yVal);
            }
            else{
                gLower.GetPoint(iPoint,yValLowUp,yValLowUp);
                if(yVal<yValLowUp){
                    gLower.SetPoint(iPoint,xVal,yVal);
                }
            }
            if(gUpper.GetN()<iPoint+1){
                gUpper.SetPoint(iPoint,xVal,yVal);
            }
            else{
                gUpper.GetPoint(iPoint,yValLowUp,yValLowUp);
                if(yVal>yValLowUp){
                    gUpper.SetPoint(iPoint,xVal,yVal);
                }
            }
            //if(ITER==DefaultIter){
            //    gDefault.SetPoint(iPoint,xVal,yVal);
            //}
            if(BestSoFarAV18){
                gBestAV18.SetPoint(iPoint,xVal,yVal);
            }
        }

        int CounterBlAV18=0;
        for(int iPoint=0; iPoint<hData->GetNbinsX(); iPoint++){
            xVal = hData->GetBinCenter(iPoint+1);
            yVal = FitBaseline->Eval(xVal);
            if(Potential==Id_AV18&&BlMin<=xVal&&xVal<=BlMax){
                if(bOuterLowerAV18.GetN()<CounterBlAV18+1){
                    bOuterLowerAV18.SetPoint(CounterBlAV18,xVal,yVal);
                }
                else{
                    bOuterLowerAV18.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal<yValLowUp){
                        bOuterLowerAV18.SetPoint(CounterBlAV18,xVal,yVal);
                    }
                }
                if(bOuterUpperAV18.GetN()<CounterBlAV18+1){
                    bOuterUpperAV18.SetPoint(CounterBlAV18,xVal,yVal);
                }
                else{
                    bOuterUpperAV18.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal>yValLowUp){
                        bOuterUpperAV18.SetPoint(CounterBlAV18,xVal,yVal);
                    }
                }
                if(BestSoFarAV18){
                    bOuterBestAV18.SetPoint(CounterBlAV18,xVal,yVal);
                }
                CounterBlAV18++;
            }
        }

        CounterBlAV18=0;
        for(int iPoint=0; iPoint<hData->GetNbinsX(); iPoint++){
            xVal = hData->GetBinCenter(iPoint+1);
            yVal = FitBaseline->Eval(xVal);
            if(Potential==Id_AV18){
                if(bLowerAV18.GetN()<CounterBlAV18+1){
                    bLowerAV18.SetPoint(CounterBlAV18,xVal,yVal);
                }
                else{
                    bLowerAV18.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal<yValLowUp){
                        bLowerAV18.SetPoint(CounterBlAV18,xVal,yVal);
                    }
                }
                if(bUpperAV18.GetN()<CounterBlAV18+1){
                    bUpperAV18.SetPoint(CounterBlAV18,xVal,yVal);
                }
                else{
                    bUpperAV18.GetPoint(iPoint,xValLowUp,yValLowUp);
                    if(yVal>yValLowUp){
                        bUpperAV18.SetPoint(CounterBlAV18,xVal,yVal);
                    }
                }
                if(BestSoFarAV18){
                    bBestAV18.SetPoint(CounterBlAV18,xVal,yVal);
                }
                CounterBlAV18++;
            }
        }
        NumAcceptedEntries++;
        delete FileGraph;

    }

    TH1D* hProjection;
    for(unsigned uPoint=0; uPoint<DevPlotNumBins; uPoint++){
        double kVal = hVarDeviationAV18->GetXaxis()->GetBinCenter(uPoint+1);
        hProjection = hVarDeviationAV18->ProjectionY("hProjection", uPoint+1,uPoint+1);
        gVarDeviationAV18_Stdv.SetPoint(uPoint,kVal,hProjection->GetMean());
        //if( (kVal>265&&kVal<295)||(uPoint==1||uPoint==2||uPoint==3||uPoint==6||uPoint==7||uPoint==8) )
        //{gVarDeviationAV18_Stdv.SetPointError(uPoint,0,20);}
        //else if( (uPoint==4||uPoint==5) )
        //{gVarDeviationAV18_Stdv.SetPointError(uPoint,0,4);}
        //else
        {gVarDeviationAV18_Stdv.SetPointError(uPoint,0,hProjection->GetStdDev());}
        delete hProjection;
    }

    printf("BestNsigma_AV18 = %.2f (chi2ndf=%.2f)\n",BestNsigma_AV18,BestChi2NdfAV18);

    printf("DefaultChi2NdfAV18=%.2f (%.2f)\n",DefaultChi2NdfAV18,DefaultNsigmaAV18);

    printf("MeanChi2NdfAV18=%.2f (%.2f)\n",MeanChi2NdfAV18,MeanNsigmaAV18);

    printf("StdvChi2NdfAV18=%.2f (%.2f)\n",StdvChi2NdfAV18,StdvNsigmaAV18);

    printf("UniformStdvChi2NdfAV18=%.2f (%.2f)\n",UniformStdvChi2NdfAV18,UniformStdvNsigmaAV18);

    printf("MinChi2NdfAV18=%.2f (%.2f)\n",MinChi2NdfAV18,MinNsigmaAV18);

    printf("MaxChi2NdfAV18=%.2f (%.2f)\n",MaxChi2NdfAV18,MaxNsigmaAV18);

    printf("--- CENTRAL INTERVAL (68%%) ---\n");
    printf(" Chi2NdfAV18: %.2f, %.2f, %.2f\n",LowCi_Chi2NdfAV18,Median_Chi2NdfAV18,UpCi_Chi2NdfAV18);
    printf(" NsigmaAV18: %.2f, %.2f, %.2f\n",LowCi_NsigmaAV18,Median_NsigmaAV18,UpCi_NsigmaAV18);

    TFile* fPlot = new TFile(FitSystFolder+"PLOT/fPlot_"+WhichBaseline+".root","recreate");
    gLowerAV18.Write();
    gUpperAV18.Write();
    gBestAV18.Write();
    hChi2NdfAV18->Write();
    hNsigmaAV18->Write();
    hVarDeviationAV18->Write();
    gVarDeviationAV18_Stdv.Write();
    hVarDevRatioAV18->Write();

    gStyle->SetCanvasPreferGL(1);
    SetStyle();
    const float right = 0.025;
    const float top = 0.025;

    TGraphErrors *grFemto_AV18;
    grFemto_AV18 = FemtoModelFitBandsSimple(&gLowerAV18, &gUpperAV18);
    grFemto_AV18->SetFillColorAlpha(fColors[2],0.3);
    grFemto_AV18->SetLineColor(fColors[2]);
    grFemto_AV18->SetLineWidth(5);
    TGraphErrors *grOuterBl_AV18;
    grOuterBl_AV18 = FemtoModelFitBandsSimple(&bOuterLowerAV18, &bOuterUpperAV18);
    grOuterBl_AV18->SetFillColorAlpha(fColors[2],0.3);
    grOuterBl_AV18->SetLineColor(fColors[2]);
    grOuterBl_AV18->SetLineWidth(5);
    gBestAV18.SetLineColor(fColors[2]);
    gBestAV18.SetLineWidth(5);
    bOuterBestAV18.SetLineColor(fColors[2]);
    bOuterBestAV18.SetLineWidth(5);
    bBestAV18.SetLineColor(fColors[2]);
    bBestAV18.SetLineWidth(4);
    bBestAV18.SetLineStyle(2);


    DLM_SubPads DlmPad(720,720);
    DlmPad.AddSubPad(0,1,0.33,1);
    DlmPad.AddSubPad(0,1,0,0.33);
    DlmPad.SetMargin(0,0.12,0.02,0.0,0.02);
    DlmPad.SetMargin(1,0.12,0.02,0.09,0.0);
    DlmPad.cd(0);

    //TCanvas *Can_CF_pL = new TCanvas("pL","pL", 0,0,720,720);
    //Can_CF_pL->SetRightMargin(right);
    //Can_CF_pL->SetTopMargin(top);

    hData->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hData->GetXaxis()->SetRangeUser(0, 350);
    hData->GetXaxis()->SetNdivisions(505);
    hData->GetYaxis()->SetRangeUser(0.8, 3.5);
    hData->SetFillColor(fFillColors[0]);
    SetStyleHisto2(hData,2,0);
    //hData->GetYaxis()->SetTitleOffset(1.0);
    hData->Draw();

    TGraphErrors *Tgraph_syserror = DrawSystematicError_FAST(hData, hSyst, fSyst, 3);
    Tgraph_syserror->SetLineColor(kWhite);

    //baselineLL->Draw("same");

    if(grFemto_AV18) {grFemto_AV18->Draw("3 same");}
    if(grOuterBl_AV18&&MODE==1) {grOuterBl_AV18->Draw("3 same");}
    if(grFemto_AV18) gBestAV18.Draw("l same");
    if(grFemto_AV18&&MODE==1) {bOuterBestAV18.Draw("l same");}
    if(grFemto_AV18&&MODE==0) {bBestAV18.Draw("l same");}
    hData->Draw("same");

    Tgraph_syserror->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror->Draw("2 same");
    //hData->Draw("pe same");

    //TString LegendSource_line1 = "Gaussian core with m_{T} scaling";
    //TString LegendSource_line1 = "Levy core with m_{T} scaling";
    //TString LegendSource_line1 = "Gaussian source";
    TString LegendSource_line1 = TString::Format("Gaussian core - fit");
    TString LegendSource_line2;
    if(WhichBaseline==0) LegendSource_line2 =  TString::Format("Constant baseline - n#sigma#in(%.1f,%.1f)",LowCi_NsigmaAV18,UpCi_NsigmaAV18);
    else if(WhichBaseline%10==1) LegendSource_line2 =  TString::Format("Linear baseline - n#sigma#in(%.1f,%.1f)",LowCi_NsigmaAV18,UpCi_NsigmaAV18);
    else if(WhichBaseline%10==2) LegendSource_line2 =  TString::Format("Quadratic baseline - n#sigma#in(%.1f,%.1f)",LowCi_NsigmaAV18,UpCi_NsigmaAV18);

    unsigned NumRows=4;
    TLegend *legend = new TLegend(0.49,0.72-0.05*NumRows,0.73,0.72);//lbrt
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(gStyle->GetTextSize()*0.90);
    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hData->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    legend->AddEntry(hCk_Fake, "p#minusp #oplus #bar{p}#minus#bar{p} pairs", "fpe");
//  legend->AddEntry(hist_CF_LL_ALAL_exp[2], "with Syst. uncertainties", "");
    //legend->AddEntry(baselineLL,"Baseline","l");
    legend->AddEntry(grFemto_AV18,"Femtoscopic fit (AV18)","l");
    legend->Draw("same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.55, 0.875, "ALICE Preliminary");
    //if(WhichDataSet==0) BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 13 TeV"));
    //else if(WhichDataSet==1) BeamText.DrawLatex(0.55, 0.825, Form("p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    //else BeamText.DrawLatex(0.55, 0.825, Form("pp #sqrt{#it{s}} = 7 TeV"));
    BeamText.DrawLatex(0.50, 0.915, "ALICE Preliminary");
    //if(DataSample=="pp13TeV_MB_Run2paper") BeamText.DrawLatex(0.50, 0.86, "pp #sqrt{#it{s}} = 13 TeV");
    //else if(DataSample.Contains("pPb")&&DataSample.Contains("5TeV")) BeamText.DrawLatex(0.50, 0.86, "p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    //else if(DataSample=="pp13TeV_HM_March19") BeamText.DrawLatex(0.50, 0.86, "pp (HM) #sqrt{#it{s}} = 13 TeV");
    //else BeamText.DrawLatex(0.50, 0.86, "ALICE pp #sqrt{#it{s}} = 7 TeV");
    BeamText.DrawLatex(0.50, 0.860, "pp (HM) #sqrt{#it{s}} = 13 TeV");

    TLatex BeamTextSource;
    BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamTextSource.SetNDC(kTRUE);
    BeamTextSource.DrawLatex(0.50, 0.805, LegendSource_line1);
    BeamTextSource.DrawLatex(0.50, 0.750, LegendSource_line2);

//INLET -------------------------------------------------------------------------------------------------------------------

    TH1F* DataHisto_Inlet = (TH1F*)hData->Clone("DataHisto_Inlet");
    DataHisto_Inlet->SetMarkerSize(hData->GetMarkerSize()*0.67);
    DataHisto_Inlet->SetLineWidth(hData->GetLineWidth()*0.67);
    DataHisto_Inlet->GetXaxis()->SetTitleSize(hData->GetXaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetLabelSize(hData->GetXaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetRangeUser(75, MODE==1?600:600);
    DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

    DataHisto_Inlet->GetYaxis()->SetTitleSize(hData->GetYaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetLabelSize(hData->GetYaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetTitleOffset(hData->GetYaxis()->GetTitleOffset()*0.67);
    DataHisto_Inlet->GetYaxis()->SetRangeUser(0.94, 1.06);

    TGraph* grFemto_AV18_Inlet = (TGraph*)grFemto_AV18->Clone("grFemto_AV18_Inlet");
    grFemto_AV18_Inlet->SetLineWidth(grFemto_AV18->GetLineWidth()*0.67);

    TGraph* gBestAV18_Inlet = (TGraph*)gBestAV18.Clone("gBestAV18_Inlet");
    gBestAV18_Inlet->SetLineWidth(gBestAV18.GetLineWidth()*0.67);

    const double fXMinInlet=0.30;
    const double fYMinInlet=0.12;
    const double fXMaxInlet=0.95;
    const double fYMaxInlet=0.50;
    TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                             fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    DataHisto_Inlet->Draw();
    //if(grFemto_AV18) grFemto_AV18_Inlet->Draw("3 same");
    //if(grFemto_NAV18) grFemto_NAV18_Inlet->Draw("3 same");
    //if(grFemto_AV18) gBestAV18_Inlet->Draw("l same");
    //if(grFemto_NAV18) gBestNAV18_Inlet->Draw("l same");
    if(grFemto_AV18) {grFemto_AV18->Draw("3 same");}
    if(grOuterBl_AV18&&MODE==1) {grOuterBl_AV18->Draw("3 same");}
    if(grFemto_AV18) gBestAV18.Draw("l same");
    if(grFemto_AV18&&MODE==1) {bOuterBestAV18.Draw("l same");}
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror->Draw("2 same");
    if(grFemto_AV18&&MODE==0) {bBestAV18.Draw("l same");}

    DlmPad.cd(1);
    TH1F* hAxis = new TH1F("hAxis", "hAxis", 35, 0, 350);
    hAxis->SetStats(false);
    hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");
    hAxis->GetXaxis()->SetRangeUser(0, 350);
    hAxis->GetYaxis()->SetRangeUser(-7.5, 7.5);
    //hData->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    //hData->GetXaxis()->SetRangeUser(0, 320);
    hAxis->GetXaxis()->SetNdivisions(505);
    //hData->GetYaxis()->SetRangeUser(0.85, 2.3);
    //hData->SetFillColor(fFillColors[0]);
    SetStyleHisto2(hAxis,2,0,2);
    //hData->GetYaxis()->SetTitleOffset(1.0);
    hAxis->Draw("");

    gVarDeviationAV18_Stdv.SetFillColorAlpha(kBlue+2,0.75);
    gVarDeviationAV18_Stdv.Draw("3,same");

    DlmPad.GetCanvas()->SaveAs(FitSystFolder+TString::Format("PLOT/Can_CF_pL_%i.pdf",WhichBaseline));

    delete hVarDeviationAV18;
    delete hVarDevRatioAV18;
    delete legend;
    delete hCk_Fake;
    delete Tgraph_syserror;
    delete grFemto_AV18;
    delete grOuterBl_AV18;
    //delete Can_CF_pL;
    delete hChi2NdfAV18;
    delete hNsigmaAV18;
    delete fPlot;
    delete DataFile;
    delete SystFile;
    delete ntFile;


}

void mT_Plots(const TString& DataSample, const bool& LevySource){

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/mT_Plots/";

    TString FileName_mT_pp;
    TString FileName_mT_pL_LO;
    TString FileName_mT_pL_NLO;
    TString FileName_mT_pXim_HAL;

    if(DataSample=="pp13TeV_HM_March19"){
        FileName_mT_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/"
        "Output/pp_FitDiff_CompareToTotal/250319/HM/Output_pp13TeV_HM_March19_McGauss_Reso_Norm.root";
        FileName_mT_pL_LO = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/"
        "Output/pp_FitDiff_CompareToTotal/250319/HM/pL/LO/Output_pp13TeV_HM_March19_McGauss_Reso_Norm.root";
        FileName_mT_pL_NLO = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/"
        "Output/pp_FitDiff_CompareToTotal/250319/HM/pL/NLO/Output_pp13TeV_HM_March19_McGauss_Reso_Norm.root";
        FileName_mT_pXim_HAL = "NOT SET";
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        FileName_mT_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/"
        "Output/pp_FitDiff_CompareToTotal/250319/pPb/Output_pPb5TeV_CPR_Mar19_McGauss_Reso_Norm.root";
        FileName_mT_pL_LO = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/"
        "Output/pp_FitDiff_CompareToTotal/250319/pPb/pL/LO/Output_pPb5TeV_CPR_Mar19_McGauss_Reso_Norm.root";
        FileName_mT_pL_NLO = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/"
        "Output/pp_FitDiff_CompareToTotal/250319/pPb/pL/NLO/Output_pPb5TeV_CPR_Mar19_McGauss_Reso_Norm.root";
        FileName_mT_pXim_HAL = "NOT SET";
    }
    else{
        printf("Trouble in Plot_mT\n");
        return;
    }

    const TString GraphName_pp = "Radius_pp";
    const TString GraphName_pL_LO = "Radius_pL";
    const TString GraphName_pL_NLO = "Radius_pL";
    const TString GraphName_pXim_HAL = "whatever";

    const TString CkBaseName_pp = "hCk";
    const TString CkBaseName_pL_LO = "hCk";
    const TString CkBaseName_pL_NLO = "hCk";
    const TString CkBaseName_pXim_HAL = "hCk";

    const TString FitBaseName_pp = "ppFit_mT";
    const TString FitBaseName_pL_LO = "pLambdaFit_mT";
    const TString FitBaseName_pL_NLO = "pLambdaFit_mT";
    const TString FitBaseName_pXim_HAL = "pXimFit_mT";

    const TString CkTotName_pp = "hTotCk";
    const TString CkTotName_pL_LO = "hTotCk";
    const TString CkTotName_pL_NLO = "hTotCk";
    const TString CkTotName_pXim_HAL = "hTotCk";

    const TString FitTotName_pp = "fit_CkTot";
    const TString FitTotName_pL_LO = "fit_CkTot";
    const TString FitTotName_pL_NLO = "fit_CkTot";
    const TString FitTotName_pXim_HAL = "fit_CkTot";

    const TString Avg_mT_Name_pp = "ppFit_TOTAL";
    const TString Avg_mT_Name_pL_LO = "pLambdaFit_TOTAL";
    const TString Avg_mTName_pL_NLO = "pLambdaFit_TOTAL";
    const TString Avg_mTName_pXim_HAL = "pXimFit_TOTAL";

    TFile* File_mT_pp=NULL;
    TFile* File_mT_pL_LO=NULL;
    TFile* File_mT_pL_NLO=NULL;
    TFile* File_mT_pXim_HAL=NULL;

    File_mT_pp = new TFile(FileName_mT_pp,"read");
    if(!File_mT_pp->IsOpen()){delete File_mT_pp; File_mT_pp=NULL;}
    printf("File_mT_pp=%p\n",File_mT_pp);
    File_mT_pL_LO = new TFile(FileName_mT_pL_LO,"read");
    if(!File_mT_pL_LO->IsOpen()){delete File_mT_pL_LO; File_mT_pL_LO=NULL;}
    printf("File_mT_pL_LO=%p\n",File_mT_pL_LO);
    File_mT_pL_NLO = new TFile(FileName_mT_pL_NLO,"read");
    if(!File_mT_pL_NLO->IsOpen()){delete File_mT_pL_NLO; File_mT_pL_NLO=NULL;}
    printf("File_mT_pL_NLO=%p\n",File_mT_pL_NLO);
    File_mT_pXim_HAL = new TFile(FileName_mT_pXim_HAL,"read");
    if(!File_mT_pXim_HAL->IsOpen()){delete File_mT_pXim_HAL; File_mT_pXim_HAL=NULL;}
    printf("File_mT_pXim_HAL=%p\n",File_mT_pXim_HAL);

    TGraph* mT_pp=NULL;
    TGraph* mT_pL_LO=NULL;
    TGraph* mT_pL_NLO=NULL;
    TGraph* mT_pXim=NULL;

    if(File_mT_pp){
        mT_pp = (TGraph*)File_mT_pp->Get(GraphName_pp);
    }
    if(File_mT_pL_LO){
        mT_pL_LO = (TGraph*)File_mT_pL_LO->Get(GraphName_pL_LO);
    }
    if(File_mT_pL_NLO){
        mT_pL_NLO = (TGraph*)File_mT_pL_NLO->Get(GraphName_pL_NLO);
    }
    if(File_mT_pXim_HAL){
        mT_pXim = (TGraph*)File_mT_pXim_HAL->Get(GraphName_pXim_HAL);
    }

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, 0.7, 2.6);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("<m_{T}> (GeV/#it{c}^{2})");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    if(LevySource) hAxis->GetYaxis()->SetTitle("Levy core (fm)");
    else hAxis->GetYaxis()->SetTitle("Gaussian core (fm)");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.90);

    //hAxis->GetXaxis()->SetNdivisions(506);
    hAxis->GetYaxis()->SetRangeUser(0.4, 1.45);

    TLegend* lLegend = new TLegend(0.6,0.65,0.95,0.95);//lbrt
    lLegend->SetName(TString::Format("lLegend"));
    lLegend->SetTextSize(0.05);
    if(mT_pp)lLegend->AddEntry(mT_pp,"p#minusp (AV18)");
    if(mT_pL_LO)lLegend->AddEntry(mT_pL_LO,"p#minus#Lambda (LO)");
    if(mT_pL_NLO)lLegend->AddEntry(mT_pL_NLO,"p#minus#Lambda (NLO)");
    if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus#Xi^{#minus} (HAL QCD)");

    TPaveText* PT1 = new TPaveText(0.4,0.9,0.975,0.975, "blNDC");//lbrt
    PT1->SetName("PT1");
    PT1->SetBorderSize(1);
    PT1->SetTextSize(0.04);
    PT1->SetFillColor(kWhite);
    PT1->SetTextFont(22);
    if(DataSample=="pp13TeV_HM_March19"){
        PT1->AddText("pp #sqrt{#it{s}} = 13 TeV (HM)");
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        PT1->AddText("p-Pb #sqrt{#it{s}_{NN}} 5.02 TeV");
    }

    TCanvas* cmT = new TCanvas("cmT", "cmT", 1);
    cmT->cd(0); cmT->SetCanvasSize(1920, 1080); cmT->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hAxis->Draw("axis");
    if(mT_pp)mT_pp->Draw("same,PL");
    if(mT_pL_LO)mT_pL_LO->Draw("same,PL");
    if(mT_pL_NLO)mT_pL_NLO->Draw("same,PL");
    if(mT_pXim)mT_pXim->Draw("same,PL");
    lLegend->Draw("same");


    gStyle->SetLineWidth(2.5);
	if(mT_pp)mT_pp->SetLineWidth(4);
	if(mT_pL_LO)mT_pL_LO->SetLineWidth(4);
	if(mT_pL_NLO)mT_pL_NLO->SetLineWidth(4);
	if(mT_pXim)mT_pXim->SetLineWidth(4);
	if(mT_pp)mT_pp->SetLineColor(kBlue+2);
	if(mT_pL_LO)mT_pL_LO->SetLineColor(kGreen+3);
	if(mT_pL_NLO)mT_pL_NLO->SetLineColor(kRed+1);
	if(mT_pXim)mT_pXim->SetLineColor(kOrange-1);
	if(mT_pp)mT_pp->SetMarkerStyle(20);
	if(mT_pL_LO)mT_pL_LO->SetMarkerStyle(21);
	if(mT_pL_NLO)mT_pL_NLO->SetMarkerStyle(22);
	if(mT_pXim)mT_pXim->SetMarkerStyle(23);
	if(mT_pp)mT_pp->SetMarkerSize(2);
	if(mT_pL_LO)mT_pL_LO->SetMarkerSize(2);
	if(mT_pL_NLO)mT_pL_NLO->SetMarkerSize(2);
	if(mT_pXim)mT_pXim->SetMarkerSize(2);
	if(mT_pp)mT_pp->SetMarkerColor(kBlue+2);
	if(mT_pL_LO)mT_pL_LO->SetMarkerColor(kGreen+3);
	if(mT_pL_NLO)mT_pL_NLO->SetMarkerColor(kRed+1);
	if(mT_pXim)mT_pXim->SetMarkerColor(kOrange-1);

    cmT->SaveAs(TString::Format("%smT_Plots_%s_%s.png",OutputFolder.Data(),DataSample.Data(),LevySource?"Levy":"Gauss"));

	gStyle->SetLineWidth(1);
	if(mT_pp)mT_pp->SetLineWidth(mT_pp->GetLineWidth()/2.5);
	if(mT_pL_LO)mT_pL_LO->SetLineWidth(mT_pL_LO->GetLineWidth()/2.5);
	if(mT_pL_NLO)mT_pL_NLO->SetLineWidth(mT_pL_NLO->GetLineWidth()/2.5);
	if(mT_pXim)mT_pXim->SetLineWidth(mT_pXim->GetLineWidth()/2.5);
    cmT->SaveAs(TString::Format("%smT_Plots_%s_%s.pdf",OutputFolder.Data(),DataSample.Data(),LevySource?"Levy":"Gauss"));
	gStyle->SetLineWidth(2.5);
	if(mT_pp)mT_pp->SetLineWidth(mT_pp->GetLineWidth()*2.5);
	if(mT_pL_LO)mT_pL_LO->SetLineWidth(mT_pL_LO->GetLineWidth()*2.5);
	if(mT_pL_NLO)mT_pL_NLO->SetLineWidth(mT_pL_NLO->GetLineWidth()*2.5);
	if(mT_pXim)mT_pXim->SetLineWidth(mT_pXim->GetLineWidth()*2.5);

    delete lLegend;
    delete hAxis;
    delete cmT;

    if(File_mT_pp){delete File_mT_pp; File_mT_pp=NULL;}
    if(File_mT_pL_LO){delete File_mT_pL_LO; File_mT_pL_LO=NULL;}
    if(File_mT_pL_NLO){delete File_mT_pL_NLO; File_mT_pL_NLO=NULL;}
    if(File_mT_pXim_HAL){delete File_mT_pXim_HAL; File_mT_pXim_HAL=NULL;}
}

int FEMTOBOYZ_MAIN(int narg, char** ARGS){
    //plotMorePreliminaries("pp");
    //plotMorePreliminaries("p#minusPb");
    //Friendship("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/SRB_111018/Full_7TeV_BL/","SystematicsRadiusBinning.root");
    //Friendship("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS2/OutputGlobal/SRB_051018/Zoomed/","SystematicsRadiusBinning_ZOOM.root");
    //plotLambda(0,"hNsigmaChi2Stdev");
    //plotLambda(3,"hNsigmaChi2Stdev");
    //plotLambda(0,"hGlobNsigma");
    //plotLambda(3,"hGlobNsigma");
    //plotLambda(0,"hGlobChi2Diff");
    //plotLambda(3,"hGlobChi2Diff");

    //plotCF("pp");

/*
    PlotDimiExclusion_ver1("/home/dmihaylov/Temp/Output/061118/MC_Default_BL_AllDataSets/","hFinalEP_Confidence",
                           "/home/dmihaylov/Temp/Output/061118/Chi2Map_SystId0.root","hNegativeLedni",3,3,
                           "pp7TeV,pp13TeV,pPb5TeV");
    PlotDimiExclusionBE_ver1("/home/dmihaylov/Temp/Output/061118/",
                             "/home/dmihaylov/Temp/Output/061118/MC_Default_BL_AllDataSets/BindingEnergy.root",
                             "/home/dmihaylov/Temp/Output/061118/MC_FullSyst_BL_AllDataSets/BindingEnergy.root",
                             "d0");
    PlotDimiExclusionBE_ver1("/home/dmihaylov/Temp/Output/061118/",
                             "/home/dmihaylov/Temp/Output/061118/MC_Default_BL_AllDataSets/BindingEnergy.root",
                             "/home/dmihaylov/Temp/Output/061118/MC_FullSyst_BL_AllDataSets/BindingEnergy.root",
                             "f0Inv");
*/
/*
    PlotDimiExclusion_ver1("/home/dmihaylov/Temp/Output/150319/MC_Default_BL_AllDataSets/","hFinalEP_Confidence",
                           "/home/dmihaylov/Temp/Output/150319/MC_Default_BL_AllDataSets/Chi2Map_SystId0.root","hNegativeLedni",3,3,
                           "pp7TeV,pp13TeV,pPb5TeV");
    PlotDimiExclusion_ver1("/home/dmihaylov/Temp/Output/150319/MC_FullSyst_BL_AllDataSets/","hFinalEP_Confidence",
                           "/home/dmihaylov/Temp/Output/150319/MC_FullSyst_BL_AllDataSets/Chi2Map_SystId0.root","hNegativeLedni",3,3,
                           "pp7TeV,pp13TeV,pPb5TeV");
    PlotDimiChi2Exclusion_ver1("/home/dmihaylov/Temp/Output/150319/MC_Default_BL_AllDataSets/","Chi2Map_SystId0.root","hChi2NdfMap",4,4,"pp7TeV,pp13TeV,pPb5TeV");
    PlotDimiExclusionBE_ver1("/home/dmihaylov/Temp/Output/150319/",
                             "/home/dmihaylov/Temp/Output/150319/MC_Default_BL_AllDataSets/BindingEnergy.root",
                             "/home/dmihaylov/Temp/Output/150319/MC_FullSyst_BL_AllDataSets/BindingEnergy.root",
                             "d0");
*/

    //void PlotLamLamFit(const TString FitFolder, const unsigned WhichDataSet, const unsigned NumSystIter,
    //               const TString DataFileName, const TString DataHistoName, const double parA, const double parB,
    //               const TString DataSystFileName)

//Data set 0: a=9.584e-01; b=1.247e-04
//Data set 0: a=9.531e-01; b=1.397e-04 (latest)
/*
Data set 0: a=9.531e-01; b=1.397e-04
         0.96; 12.22 (r=1.15)
  at 100.00 ck is 0.87
  pp0 = 0.96; pp1 = 12.22
Data set 1: a=9.997e-01; b=-7.879e-06
         -1000000.00; -1000000.00 (r=1.40)
  at 100.00 ck is 0.93
  pp0 = 0.96; pp1 = 12.22
Data set 2: a=1.003e+00; b=-2.650e-06
         -1000000.00; -1000000.00 (r=1.12)
  at 100.00 ck is 0.87
  pp0 = 0.96; pp1 = 12.22
*/


/*
    PlotLamLamFit("/home/dmihaylov/Temp/Output/150319/StuffForLamLamFit/",0,27,
                  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root",
                  "hCkTotNormWeight",9.531e-01,1.397e-04,
                  "/home/dmihaylov/CernBox/Femto_pp13/data_old/C2totalsysLL.root");

    PlotLamLamFit("/home/dmihaylov/Temp/Output/150319/StuffForLamLamFit/",1,27,
                  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root",
                  "hCkTotNormWeight",9.997e-01,-7.879e-06,
                  "/home/dmihaylov/CernBox/SystematicsAndCalib/pPbRun2_MB/C2totalsysLL.root");

    PlotLamLamFit("/home/dmihaylov/Temp/Output/150319/StuffForLamLamFit/",2,27,
                  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/301018/CFOutput_LL.root",
                  "hCkTotReweightMeV",1.003e+00,-2.650e-06,
                  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Andi/Systematics/C2totalsysLL.root");

    plotManyLambdaLambdaModels("/home/dmihaylov/Temp/Output/150319/StuffForLamLamFit/","pp");
    plotManyLambdaLambdaModels("/home/dmihaylov/Temp/Output/150319/StuffForLamLamFit/","pPb");
*/

    //mT_Plots("pp13TeV_HM_March19",false);
    //mT_Plots("pPb5TeV_CPR_Mar19",false);

    //Plot_pL_FASTsyst(0);
    //Plot_pL_FASTsyst(11);
    //Plot_pL_FASTsyst(12);

    Plot_pp_FASTsyst(0);
    Plot_pp_FASTsyst(1);
    Plot_pp_FASTsyst(2);
    Plot_pp_FASTsyst(11);
    Plot_pp_FASTsyst(12);

    return 0;
}
