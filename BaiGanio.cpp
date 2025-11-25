
#include "BaiGanio.h"

#include "TRandom3.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TGraphErrors.h"

#include "CATS.h"
#include "CATStools.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "EnvVars.h"
#include "CATSconstants.h"
#include "CommonAnaFunctions.h"

#include "TROOT.h"

#include "DLM_Histo.h"



#include <fstream>   // at the top of your file
#include <iomanip>

#include "Math/MinimizerOptions.h"



//4D: //f0,d0,kstar,rstar
DLM_Histo<float>* dlm_wf_u_pp_coarse = NULL;
DLM_Histo<float>* dlm_wf_u_pp_fine = NULL;

struct Ganio_pp_pot_pars
{
    float V0,mu0,V1,mu1,f0,d0;
    int num_entries=0;
    //void add_entry();
  Ganio_pp_pot_pars(){
    Ganio_pp_pot_pars(0);
  }
   Ganio_pp_pot_pars(double value){
    *this = 0;
    }
  

  bool operator+=(const Ganio_pp_pot_pars& other){
    //double w = double(num_entries)/double(num_entries+other.num_entries);
    //V0 = V0*w + other.V0*(1.-w);
    //mu0 = mu0*w + other.mu0*(1.-w);
    //V1 = V1*w + other.V1*(1.-w);
    //mu1 = mu1*w + other.mu1*(1.-w);
    //f0 = f0*w + other.f0*(1.-w);
    //d0 = d0*w + other.d0*(1.-w);
    //num_entries += other.num_entries;

    V0 += other.V0;
    mu0 += other.mu0;
    V1 += other.V1;
    mu1 += other.mu1;
    f0 += other.f0;
    d0 += other.d0;
    num_entries += other.num_entries;

    return true;
  }
  bool operator/=(const double& value){
    V0 /= value;
    mu0 /= value;
    V1 /= value;
    mu1 /= value;
    f0 /= value;
    d0 /= value;
    num_entries = TMath::Nint(double(num_entries)/value);
    return true;
  }

    bool operator=(const Ganio_pp_pot_pars& other){
        V0 = other.V0;
        mu0 = other.mu0;
        V1 = other.V1;
        mu1 = other.mu1;
        f0 = other.f0;
        d0 = other.d0;
        num_entries = other.num_entries;
        return true;
  }
    bool operator=(const double& value){
        V0 = value;
        mu0 = value;
        V1 = value;
        mu1 = value;
        f0 = value;
        d0 = value;
        num_entries = 0;
        return true;
  }
  Ganio_pp_pot_pars operator*(const double& value){
      Ganio_pp_pot_pars Result;
        Result.V0 = V0 * value;
        Result.mu0 = mu0 * value;
        Result.V1 = V1 * value;
        Result.mu1 = mu1 * value;
        Result.f0 = f0 * value;
        Result.d0 = d0 * value;
        Result.num_entries = num_entries;
        return Result;
  }
  void compute_average(){
    *this /= num_entries;
  }
  void print(){
    printf("%p: f0=%.3e d0=%.3e V0=%.3e mu0=%.3e V1=%.3e mu2=%.3e ne=%i\n",this,f0,d0,V0,mu0,V1,mu1,num_entries);
  }
};


void SetUp_Ganio_pp_pot(const int n_f0, const double min_f0, const double max_f0,
                        const int n_d0, const double min_d0, const double max_d0,
                        const int n_kstar, const double min_kstar, const double max_kstar, 
                        const int n_rstar, const double min_rstar, const double max_rstar,
                        TNtuple& input_ntuple, const int pot_id, DLM_Histo<float>*& dlm_wf_u_pp){

    printf("Hello there, lets set things up!\n");

    if(dlm_wf_u_pp) delete dlm_wf_u_pp;
    dlm_wf_u_pp = new DLM_Histo<float> ();
    //f0,d0,kstar,rstar
    dlm_wf_u_pp->SetUp(4);
    dlm_wf_u_pp->SetUp(0, n_f0, min_f0, max_f0);
    dlm_wf_u_pp->SetUp(1, n_d0, min_d0, max_d0);
    dlm_wf_u_pp->SetUp(2, n_kstar, min_kstar, max_kstar);
    dlm_wf_u_pp->SetUp(3, n_rstar, min_rstar, max_rstar);
    dlm_wf_u_pp->Initialize();
    
    //f0,d0
    DLM_Histo<Ganio_pp_pot_pars> dlm_par_map;
    dlm_par_map.SetUp(2);
    dlm_par_map.SetUp(0, n_f0, min_f0, max_f0);
    dlm_par_map.SetUp(1, n_d0, min_d0, max_d0);   
    dlm_par_map.Initialize();
    
    float SEED, POT_ID, RedM, f0, d0, V0, mu0, V1, mu1;
    input_ntuple.SetBranchAddress("SEED",&SEED);
    input_ntuple.SetBranchAddress("POT_ID",&POT_ID);
    input_ntuple.SetBranchAddress("RedM",&RedM);
    input_ntuple.SetBranchAddress("f0",&f0);
    input_ntuple.SetBranchAddress("d0",&d0);
    input_ntuple.SetBranchAddress("V0",&V0);
    input_ntuple.SetBranchAddress("mu0",&mu0);
    input_ntuple.SetBranchAddress("V1",&V1);
    input_ntuple.SetBranchAddress("mu1",&mu1);
    
    
    gROOT->cd();
    TGraphErrors* g_visual_mu0 = new TGraphErrors();
    g_visual_mu0->SetName("g_visual_mu0");

    //I did a visual cut to get rid of offliers. Simply check if mu < some pol2 function, and hard cut on V
    TF1* fLowerLimit_mu0 = new TF1("fLowerLimit_mu0", "[0]+[1]*x+[2]*x*x",-128,0);
    TF1* fUpperLimit_mu0 = new TF1("fUpperLimit_mu0", "[0]+[1]*x+[2]*x*x",-128,0);

    TF1* fLowerLimit_V0 = new TF1("fLowerLimit_V0", "[0]+[1]*x+[2]*x*x",0,8);
    TF1* fUpperLimit_V0 = new TF1("fUpperLimit_V0", "[0]+[1]*x+[2]*x*x",0,8);

    TH2F* h_V0_mu0_input = new TH2F("h_V0_mu0_input","h_V0_mu0_input",512,-128,0,512,0,8);  
    TH2F* h_V0_mu0_selection = new TH2F("h_V0_mu0_selection","h_V0_mu0_selection",512,-128,0,512,0,8);  

    //cut out regions in the parameter space (inspect visually)
    if(pot_id==1){

        
        g_visual_mu0->SetPoint(0, -85, 1.10);
        g_visual_mu0->SetPointError(0, 0, 0.1);

        g_visual_mu0->SetPoint(1, -60, 1.29);
        g_visual_mu0->SetPointError(1, 0, 0.1);  

        g_visual_mu0->SetPoint(2, -40, 1.57);
        g_visual_mu0->SetPointError(2, 0, 0.1); 
        
        g_visual_mu0->SetPoint(3, -25, 1.96);
        g_visual_mu0->SetPointError(3, 0, 0.1);
        
        g_visual_mu0->SetPoint(4, -17, 2.34);
        g_visual_mu0->SetPointError(4, 0, 0.1);   
        
        g_visual_mu0->Fit(fLowerLimit_mu0,"S, N, R, M");

        fLowerLimit_mu0->SetParameter(0, fLowerLimit_mu0->GetParameter(0)-0.2);

        fUpperLimit_mu0->SetParameter(0, fLowerLimit_mu0->GetParameter(0) + 0.2*2);
        fUpperLimit_mu0->SetParameter(1, fLowerLimit_mu0->GetParameter(1));
        fUpperLimit_mu0->SetParameter(2, fLowerLimit_mu0->GetParameter(2));

        //fLowerLimit_mu0->SetParameter(0, 3.025);
        //fLowerLimit_mu0->SetParameter(1, 0.05434);
        //fLowerLimit_mu0->SetParameter(2, 0.0003432);

        fLowerLimit_V0->SetParameter(0, -88);
        fLowerLimit_V0->SetParameter(1, 0);
        fLowerLimit_V0->SetParameter(2, 0);        

        //fUpperLimit_mu0->SetParameter(0, 3.325);
        //fUpperLimit_mu0->SetParameter(1, 0.05434);
        //fUpperLimit_mu0->SetParameter(2, 0.0003432);

        fUpperLimit_V0->SetParameter(0, -16);
        fUpperLimit_V0->SetParameter(1, 0);
        fUpperLimit_V0->SetParameter(2, 0);        
    }

    if(pot_id==11){

        
        g_visual_mu0->SetPoint(0, -39, 1.55);
        g_visual_mu0->SetPointError(0, 0, 0.1);

        g_visual_mu0->SetPoint(1, -32, 1.7);
        g_visual_mu0->SetPointError(1, 0, 0.1);  

        g_visual_mu0->SetPoint(2, -27, 1.85);
        g_visual_mu0->SetPointError(2, 0, 0.1); 
        
        g_visual_mu0->SetPoint(3, -20.75, 2.1);
        g_visual_mu0->SetPointError(3, 0, 0.1);
        
        g_visual_mu0->SetPoint(4, -15.5, 2.4);
        g_visual_mu0->SetPointError(4, 0, 0.1);   

        g_visual_mu0->SetPoint(5, -10, 2.95);
        g_visual_mu0->SetPointError(5, 0, 0.1);   

        g_visual_mu0->SetPoint(6, -7.5, 3.35);
        g_visual_mu0->SetPointError(6, 0, 0.1);   


        g_visual_mu0->Fit(fLowerLimit_mu0,"S, N, R, M");

        fLowerLimit_mu0->SetParameter(0, fLowerLimit_mu0->GetParameter(0)-0.2);

        fUpperLimit_mu0->SetParameter(0, fLowerLimit_mu0->GetParameter(0) + 0.2*2);
        fUpperLimit_mu0->SetParameter(1, fLowerLimit_mu0->GetParameter(1));
        fUpperLimit_mu0->SetParameter(2, fLowerLimit_mu0->GetParameter(2));

        //fLowerLimit_mu0->SetParameter(0, 3.025);
        //fLowerLimit_mu0->SetParameter(1, 0.05434);
        //fLowerLimit_mu0->SetParameter(2, 0.0003432);

        fLowerLimit_V0->SetParameter(0, -26);
        fLowerLimit_V0->SetParameter(1, 0);
        fLowerLimit_V0->SetParameter(2, 0);        

        //fUpperLimit_mu0->SetParameter(0, 3.325);
        //fUpperLimit_mu0->SetParameter(1, 0.05434);
        //fUpperLimit_mu0->SetParameter(2, 0.0003432);

        fUpperLimit_V0->SetParameter(0, -9);
        fUpperLimit_V0->SetParameter(1, 0);
        fUpperLimit_V0->SetParameter(2, 0);        
    }
    if(pot_id==201){
       printf("pot_id==201 not set up yet !!!!\n");
    }
    if(pot_id==202){
       printf("pot_id==202 not set up yet !!!!\n");
    }
    if(pot_id==203){
       printf("pot_id==203 not set up yet !!!!\n");
    }
    if(pot_id==204){
       printf("pot_id==204 not set up yet !!!!\n");
    }
    if(pot_id==205){
       printf("pot_id==205 not set up yet !!!!\n");
    }

    for(unsigned uEntry=0; uEntry<input_ntuple.GetEntries(); uEntry++){
        if(uEntry%1000000==0) printf("%u of %u M\n", uEntry/1000000, input_ntuple.GetEntries()/1000000);
        input_ntuple.GetEntry(uEntry);
        if(TMath::Nint(POT_ID)!=pot_id) continue;

        if(f0<min_f0 || f0>max_f0) continue;
        if(d0<min_d0 || d0>max_d0) continue;

        h_V0_mu0_input->Fill(V0,mu0);

        if(V0<fLowerLimit_V0->Eval(V0)) continue;
        if(V0>fUpperLimit_V0->Eval(V0)) continue;
        if(mu0<fLowerLimit_mu0->Eval(V0)) continue;
        if(mu0>fUpperLimit_mu0->Eval(V0)) continue;

        h_V0_mu0_selection->Fill(V0,mu0);

        //превъртяхме потенциала
        //if(pot_id==11&&V0<-32) continue;

        //THINK HOW TO DO THE APPROXIMATE ESTIMATION OF THE V0/mu0 <-> f0,d0, WE ACTUALLY JUST NEED A SINGLE SOLUTION PER f0/d0 pair
        //use dlm_par_map to get some averages
        Ganio_pp_pot_pars ParCombo;
        
        ParCombo.V0 = V0;
        ParCombo.mu0 = mu0;
        ParCombo.V1 = V1;
        ParCombo.mu1 = mu1;
        ParCombo.f0 = f0;
        ParCombo.d0 = d0;
        ParCombo.num_entries = 1;

        dlm_par_map.AddAt(f0,d0,ParCombo);
    }

    //printf("WTF %u\n",dlm_par_map.GetNbins());
    for(unsigned uBin=0; uBin<dlm_par_map.GetNbins(); uBin++){
        //Ganio_pp_pot_pars& element = dlm_par_map.GetBinElement(uBin);
        //element.compute_average();
        //dlm_par_map.GetBinElement(uBin).print();
        dlm_par_map.GetBinElement(uBin).compute_average();
        //dlm_par_map.GetBinElement(uBin).print();
        unsigned WhichBin[2];
        dlm_par_map.GetBinCoordinates(uBin, WhichBin);
        //printf("  %f %f\n",dlm_par_map.GetBinLowEdge(0,WhichBin[0]),dlm_par_map.GetBinUpEdge(0,WhichBin[0]));
        //printf("   %f %f\n",dlm_par_map.GetBinLowEdge(1,WhichBin[1]),dlm_par_map.GetBinUpEdge(1,WhichBin[1]));
        if( dlm_par_map.GetBinLowEdge(0,WhichBin[0])<=17.3 && dlm_par_map.GetBinUpEdge(0,WhichBin[0])>17.3 && 
            dlm_par_map.GetBinLowEdge(1,WhichBin[1])<=2.8 && dlm_par_map.GetBinUpEdge(1,WhichBin[1])>2.8){

            dlm_par_map.GetBinElement(uBin).print();
            //printf(" bc %.3ef %.3ef\n",dlm_par_map.GetBinCenter(0,WhichBin[0]),dlm_par_map.GetBinCenter(1,WhichBin[1]));
            //printf("\n");
            //usleep(40e3);
        }

    }

    printf("The parameters for pp:\n");
    dlm_par_map.Eval2D(17.3,2.8,false).print();

    const double kMin = 0;
    const double kMax = 60;
    const unsigned kSteps = 10;

    const double SourceSize = 1.0;

    printf("xCheck:\n");
    CATS Kitty;
    Kitty.SetMomBins(kSteps, kMin, kMax);

    CATSparameters sPars(CATSparameters::tSource,1,true);
    Kitty.SetAnaSource(GaussSource, sPars);
    Kitty.SetAnaSource(0, SourceSize);
    Kitty.SetUseAnalyticSource(true);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetNormalizedSource(false);
    Kitty.SetAutoNormSource(true);

    Kitty.SetRedMass(0.5*Mass_p);
    Kitty.SetQ1Q2(0);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetNumChannels(1);
    Kitty.SetSpin(0,0);
    Kitty.SetNumPW(0,1);
    Kitty.SetChannelWeight(0, 1);

    Kitty.SetEpsilonConv(2e-9);
    Kitty.SetEpsilonProp(2e-9);
    
    
    
    CATS FatCat;
    FatCat.SetMomBins(n_kstar, min_kstar, max_kstar);
    FatCat.SetAnaSource(GaussSource, sPars);
    FatCat.SetAnaSource(0, SourceSize);
    FatCat.SetUseAnalyticSource(true);

    FatCat.SetMomentumDependentSource(false);
    FatCat.SetThetaDependentSource(false);
    FatCat.SetNormalizedSource(false);
    FatCat.SetAutoNormSource(true);

    FatCat.SetRedMass(0.5*Mass_p);
    FatCat.SetQ1Q2(1);
    FatCat.SetQuantumStatistics(true);
    FatCat.SetNumChannels(2);
    FatCat.SetSpin(0,0);
    FatCat.SetSpin(1,1);
    FatCat.SetNumPW(0,1);
    FatCat.SetChannelWeight(0, 0.25);
    FatCat.SetChannelWeight(1, 0.75);

    FatCat.SetNotifications(CATS::nWarning);

    FatCat.SetEpsilonConv(2e-9);
    FatCat.SetEpsilonProp(2e-9);
    //FatCat.SetMaxRad(128);



    if(pot_id==11){
        double PotPars[2]={dlm_par_map.Eval2D(17.3,2.8,false).V0, dlm_par_map.Eval2D(17.3,2.8,false).mu0};
//PotPars[0] = -14.08;
//PotPars[1] = 2.547;
        CATSparameters cPotPars(CATSparameters::tPotential,2,true);
        cPotPars.SetParameters(PotPars);
        Kitty.SetShortRangePotential(0,0,SquareWell,cPotPars);
        FatCat.SetShortRangePotential(0,0,SquareWell,cPotPars);
        printf(" sw pot_pars: %.3f %3f\n", Kitty.GetPotPar(0,0,0), Kitty.GetPotPar(0,0,1));
    }
    else if(pot_id==1){
        double PotPars[2]={dlm_par_map.Eval2D(17.3,2.8,false).V0, dlm_par_map.Eval2D(17.3,2.8,false).mu0};
        CATSparameters cPotPars(CATSparameters::tPotential,2,true);
        cPotPars.SetParameters(PotPars);
        Kitty.SetShortRangePotential(0,0,SingleGauss,cPotPars);
        FatCat.SetShortRangePotential(0,0,SingleGauss,cPotPars);
        printf(" sg pot_pars: %.3f %3f\n", Kitty.GetPotPar(0,0,0), Kitty.GetPotPar(0,0,1));        
    }
    //DG
    else if(pot_id>=200 && pot_id<=299){
        double PotPars[4]={dlm_par_map.Eval2D(17.3,2.8,false).V0, dlm_par_map.Eval2D(17.3,2.8,false).mu0,
                            dlm_par_map.Eval2D(17.3,2.8,false).V1, dlm_par_map.Eval2D(17.3,2.8,false).mu1};
        CATSparameters cPotPars(CATSparameters::tPotential,4,true);
        cPotPars.SetParameters(PotPars);
        Kitty.SetShortRangePotential(0,0,DoubleGaussSum,cPotPars);
        FatCat.SetShortRangePotential(0,0,DoubleGaussSum,cPotPars);
        printf(" dg pot_pars: %.3f %3f\n", Kitty.GetPotPar(0,0,0), Kitty.GetPotPar(0,0,1), Kitty.GetPotPar(0,0,2), Kitty.GetPotPar(0,0,3));     
    }

    printf("we have a cat\n");
    Kitty.KillTheCat(); 

    for(unsigned uf0=0; uf0<n_f0; uf0++){
        unsigned WhichBin[4];
        WhichBin[0] = uf0;
        double f0_val = dlm_wf_u_pp->GetBinCenter(0,uf0);
        printf("uf0 %u of %u (%.2f)\n",uf0,n_f0,f0_val);
        for(unsigned ud0=0; ud0<n_d0; ud0++){
            WhichBin[1] = ud0;
            double d0_val = dlm_wf_u_pp->GetBinCenter(1,ud0);
            //printf("    ud0 %u of %u (%.2f)\n",ud0,n_d0,d0_val);
            double PotPars[2]={dlm_par_map.Eval2D(f0_val,d0_val,false).V0, dlm_par_map.Eval2D(f0_val,d0_val,false).mu0};
            //printf("    V = %.3f; mu = %.4f\n",PotPars[0],PotPars[1]);
            //dlm_par_map.GetBinContent(uf0,ud0).print();
            //usleep(100e3);
            FatCat.SetShortRangePotential(0,0,0,PotPars[0]);
            FatCat.SetShortRangePotential(0,0,1,PotPars[1]);
            FatCat.KillTheCat(); 
            //put the vals for the WF
            for(unsigned uks=0; uks<n_kstar; uks++){
                WhichBin[2] = uks;
                for(unsigned urs=0; urs<n_rstar; urs++){
                    WhichBin[3] = urs;
                    double rad = dlm_wf_u_pp->GetBinCenter(3, urs);
                    //4D: //f0,d0,kstar,rstar
                    //DLM_Histo<float>* dlm_wf_u_pp = NULL;
                    //dlm_wf_u_pp->SetUp(0, n_f0, min_f0, max_f0);
                    //dlm_wf_u_pp->SetUp(1, n_d0, min_d0, max_d0);
                    //dlm_wf_u_pp->SetUp(2, n_kstar, min_kstar, max_kstar);
                    //dlm_wf_u_pp->SetUp(3, n_rstar, min_rstar, max_rstar);
                    dlm_wf_u_pp->SetBinContent(WhichBin,FatCat.EvalRadialWaveFunction(uks, 0, 0, rad, false).real()*hbarc);
                    //printf("%u %u %u %u %f vs %f\n", 
                    //    WhichBin[0],WhichBin[1],WhichBin[2],WhichBin[3], FatCat.EvalRadialWaveFunction(uks, 0, 0, rad, false).real(),
                    //    WhichBin[0],WhichBin[1],WhichBin[2],WhichBin[3], FatCat.EvalReferenceRadialWF(uks, 0, rad, false).real());
                    //usleep(50e3);
                }
            }
        }
    }

    //dlm_wf_u_pp->SetUp(0, n_f0, min_f0, max_f0);
    //dlm_wf_u_pp->SetUp(1, n_d0, min_d0, max_d0);
    //dlm_wf_u_pp->SetUp(2, n_kstar, min_kstar, max_kstar);
    //dlm_wf_u_pp->SetUp(3, n_rstar, min_rstar, max_rstar);


    double ScatLen, EffRan;
    TH1F *hFit;
    TF1 *fitSP;
    GetScattParameters(Kitty, ScatLen, EffRan, hFit, fitSP, 2, false, false, 0);    
    printf("I got from the cat f0 = %.2f and d0 = %.2f\n", ScatLen, EffRan);
    printf("ERE set up (q1q2): %f\n", fitSP->GetParameter(0));

    printf("PARS: %f %f\n", Kitty.GetPotPar(0,0,0), Kitty.GetPotPar(0,0,1));
    FatCat.SetShortRangePotential(0,0,0,Kitty.GetPotPar(0,0,0));
    FatCat.SetShortRangePotential(0,0,1,Kitty.GetPotPar(0,0,1));
    FatCat.KillTheCat(CATS::kAllChanged);

    DLM_CommonAnaFunctions AnalysisObject; 
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

    CATS Kitty_AV18;
    Kitty_AV18.SetMomBins(n_kstar, min_kstar, max_kstar);
    Kitty_AV18.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pp(Kitty_AV18,"AV18_s","Gauss",0,0);//McLevyNolan_Reso
    Kitty_AV18.SetAnaSource(0,FatCat.GetAnaSourcePar(0));
    Kitty_AV18.KillTheCat();

    TFile fOut(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/GanioDump.root", GetFemtoOutputFolder()),"recreate");

    TGraph gKitty;
    gKitty.SetName("gKitty");
    TGraph gFatCat;
    gFatCat.SetName("gFatCat");
    TGraph gAv18;
    gAv18.SetName("gAv18");

    TGraph gFatCatPS;
    gFatCatPS.SetName("gFatCatPS");
    TGraph gAv18PS;
    gAv18PS.SetName("gAv18PS");

    TGraph gFatCatU50;
    gFatCatU50.SetName("gFatCatU50");
    TGraph gAv18U50;
    gAv18U50.SetName("gAv18U50");    

    TGraphErrors gPotAv18;
    gPotAv18.SetName("gPotAv18");
    for(double frad=0.05; frad<10; frad+=0.05){
        gPotAv18.SetPoint(gPotAv18.GetN(), frad, Kitty_AV18.EvaluateThePotential(0,0,10,frad));
        gPotAv18.SetPointError(gPotAv18.GetN()-1, 0.0, 0.1);
    }

    for(unsigned uKstar=0; uKstar<kSteps; uKstar++){
        double kstar, Ck_Kitty;
        kstar = Kitty.GetMomentum(uKstar);
        Ck_Kitty = Kitty.GetCorrFun(uKstar);
        gKitty.SetPoint(uKstar, kstar, Ck_Kitty); 
    }
    for(unsigned uKstar=0; uKstar<n_kstar; uKstar++){
        double kstar, Ck_FatCat, Ck_Av18;
        kstar = FatCat.GetMomentum(uKstar);
        Ck_FatCat = FatCat.GetCorrFun(uKstar);
        Ck_Av18 = Kitty_AV18.GetCorrFun(uKstar);
        gFatCat.SetPoint(uKstar, kstar, Ck_FatCat);
        gAv18.SetPoint(uKstar, kstar, Ck_Av18);

        gFatCatPS.SetPoint(uKstar, kstar, FatCat.GetPhaseShift(uKstar,0,0));
        gAv18PS.SetPoint(uKstar, kstar, Kitty_AV18.GetPhaseShift(uKstar,0,0));        
    }
    for(double rad=0.025; rad<=8; rad+=0.05){
        gFatCatU50.SetPoint(gFatCatU50.GetN(), rad, FatCat.EvalRadialWaveFunction(FatCat.GetMomBin(50), 0, 0, rad, false).real());
        gAv18U50.SetPoint(gAv18U50.GetN(), rad, Kitty_AV18.EvalRadialWaveFunction(Kitty_AV18.GetMomBin(50), 0, 0, rad, false).real());
    }

    gKitty.Write();
    gFatCat.Write();
    gAv18.Write();

    gFatCatPS.Write();
    gAv18PS.Write();    
    gPotAv18.Write();

    gFatCatU50.Write();
    gAv18U50.Write();

    //fitSP->SetParameter(1, 17.3);
    //fitSP->SetParameter(2, 2.8);
    hFit->Write();
    fitSP->Write();

    fLowerLimit_mu0->Write();
    fUpperLimit_mu0->Write();

    fLowerLimit_V0->Write();
    fUpperLimit_V0->Write();

    h_V0_mu0_input->Write();
    h_V0_mu0_selection->Write();

    g_visual_mu0->Write();

    delete fLowerLimit_mu0;
    delete fUpperLimit_mu0;
    delete fLowerLimit_V0;
    delete fUpperLimit_V0;

    delete h_V0_mu0_input;
    delete h_V0_mu0_selection;
    delete g_visual_mu0;
}

double Ganio_pp_pot(double* pars){

}


//NumEvents in K
void Generate_GaussPotPars_For_Kali(bool ASCII, TString Description, TString PotType, int PotVar, float MinRedM, float MaxRedM, unsigned NumEvents, int SEED){

    NumEvents *= 1000;

    //const double MinRedM = 450;
    //const double MaxRedM = 600;
    
    unsigned* NumBinsPar = new unsigned [4];
    float* MinPar = new float [4];
    float* MaxPar = new float [4];
    bool* LogSampling = new bool [4];

    //double MinAmp1;
    //double MaxAmp1;
    //unsigned NumAmp1;

    //double MinAmp2;
    //double MaxAmp2;
    //unsigned NumAmp2;

    //double MinRng1;
    //double MaxRng1;
    //unsigned NumRng1;

    //double MinRng2 = 1.0;
    //double MaxRng2 = 2.1;
    //unsigned NumRng2;

    const double kMin = 0;
    const double kMax = 60;
    const unsigned kSteps = 6;

    //17.3
    double f0_goal_min = 17.3*0.7;
    double f0_goal_max = 17.3*1.3;
    unsigned Num_f0 = 32;
//f0_goal_min = 17.3*0.85;//temp
//f0_goal_max = 17.3*1.15;//temp
//f0_goal_min = 17.0;
//f0_goal_max = 20.0;
//Num_f0 = 8;
    //2.8
    double d0_goal_min = 2.8*0.7;
    double d0_goal_max = 2.8*1.3;
    unsigned Num_d0 = 32;//was 16
//d0_goal_min = 2.8*0.85;//temp
//d0_goal_max = 2.8*1.15;//temp
//d0_goal_min = 2.5;
//d0_goal_max = 2.7;
//Num_d0 = 8;
    

    unsigned EachXevent_rnd = 4;

    //this turned out VERY tricky, we need fine, really fine binning
    if(PotType=="SquareWell"){
        Num_f0 = 48;
        Num_d0 = 128;
        //the values here are too wild, so less random
        EachXevent_rnd = 8;
        //f0_goal_min = 17.3*0.8;//temp
        //f0_goal_max = 17.3*1.4;//temp
    }

    TRandom3 rangen(SEED);

    CATS Kitty;
    Kitty.SetMomBins(kSteps, kMin, kMax);

    CATSparameters sPars(CATSparameters::tSource,1,true);
    Kitty.SetAnaSource(GaussSource, sPars);
    Kitty.SetAnaSource(0, 1.0);
    Kitty.SetUseAnalyticSource(true);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetNormalizedSource(false);
    Kitty.SetAutoNormSource(true);

    Kitty.SetRedMass(0.5*Mass_p);
    Kitty.SetQ1Q2(0);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetNumChannels(1);
    Kitty.SetSpin(0,0);
    Kitty.SetNumPW(0,1);
    Kitty.SetChannelWeight(0, 1);
    double PotPars[2]={0,1};
    CATSparameters cPotPars(CATSparameters::tPotential,4,true);
    cPotPars.SetParameters(PotPars);

    printf("we have a cat\n");



    float potential_id = 0;
    unsigned NumPotPars = 0;

    double epsc = 1e-8;
    double epsp = 1e-8;
    
    if(PotType=="Gauss"){
        potential_id = 1;
        NumPotPars = 2;

        Kitty.SetShortRangePotential(0,0,SingleGauss,cPotPars);
        MinPar[0] = -90;
        MaxPar[0] = -12;
        NumBinsPar[0] = 35;
        LogSampling[0] = false;

        MinPar[1] = 1.0;
        MaxPar[1] = 2.6;
        NumBinsPar[1] = 44;
        LogSampling[1] = false;
    }
    else if(PotType=="DoubleGauss"){
//do these variations for the core:
//3000, 0.4
//4500, 0.4
//2250, 0.4
//4500, 0.3
//2250, 0.6

        potential_id = 200;
        potential_id += PotVar;
        NumPotPars = 4;

        Kitty.SetShortRangePotential(0,0,DoubleGaussSum,cPotPars);

        if(PotVar==1){

            MinPar[0] = -420;//480
            MaxPar[0] = -30;//160
            NumBinsPar[0] = 24*4;
            LogSampling[0] = false;
            MinPar[1] = 0.8;//0.7
            MaxPar[1] = 2.0;//1.4
            NumBinsPar[1] = 24*4;
            LogSampling[1] = false;

            MinPar[2] = 3000;
            MaxPar[2] = 3000;
            NumBinsPar[2] = 1;
            LogSampling[2] = true;
            MinPar[3] = 0.4;
            MaxPar[3] = 0.4;
            NumBinsPar[3] = 1;
            LogSampling[3] = false;            
        }
        else if(PotVar==2){

            MinPar[0] = -500;//560
            MaxPar[0] = -60;//240
            NumBinsPar[0] = 24*4;
            LogSampling[0] = false;
            MinPar[1] = 0.75;//0.6
            MaxPar[1] = 1.8;//1.3
            NumBinsPar[1] = 24*4;
            LogSampling[1] = false;

            MinPar[2] = 4500;
            MaxPar[2] = 4500;
            NumBinsPar[2] = 1;
            LogSampling[2] = true;
            MinPar[3] = 0.4;
            MaxPar[3] = 0.4;
            NumBinsPar[3] = 1;
            LogSampling[3] = false;            
        }
        else if(PotVar==3){


            MinPar[0] = -300;//420
            MaxPar[0] = -35;//120
            NumBinsPar[0] = 24*4;
            LogSampling[0] = false;
            MinPar[1] = 0.85;//0.7
            MaxPar[1] = 1.95;//1.4
            NumBinsPar[1] = 24*4;
            LogSampling[1] = false;

            MinPar[2] = 2250;
            MaxPar[2] = 2250;
            NumBinsPar[2] = 1;
            LogSampling[2] = true;
            MinPar[3] = 0.4;
            MaxPar[3] = 0.4;
            NumBinsPar[3] = 1;
            LogSampling[3] = false;            
        }        
        else if(PotVar==4){

            MinPar[0] = -260;//480
            MaxPar[0] = -25;//160
            NumBinsPar[0] = 24*4;
            LogSampling[0] = false;
            MinPar[1] = 0.8;//0.7
            MaxPar[1] = 2.2;//1.4
            NumBinsPar[1] = 24*4;
            LogSampling[1] = false;

            MinPar[2] = 4500;
            MaxPar[2] = 4500;
            NumBinsPar[2] = 1;
            LogSampling[2] = true;
            MinPar[3] = 0.3;
            MaxPar[3] = 0.3;
            NumBinsPar[3] = 1;
            LogSampling[3] = false;            
        }  
        else if(PotVar==5){
            MinPar[0] = -1600;//480
            MaxPar[0] = -90;//160
            NumBinsPar[0] = 24*4;
            LogSampling[0] = false;
            MinPar[1] = 0.7;//0.7
            MaxPar[1] = 1.6;//1.4
            NumBinsPar[1] = 24*4;
            LogSampling[1] = false;

            MinPar[2] = 2250;
            MaxPar[2] = 2250;
            NumBinsPar[2] = 1;
            LogSampling[2] = true;
            MinPar[3] = 0.6;
            MaxPar[3] = 0.6;
            NumBinsPar[3] = 1;
            LogSampling[3] = false;            
        }
        else{

            MinPar[0] = -100;
            MaxPar[0] = 0;
            NumBinsPar[0] = 24;
            LogSampling[0] = false;
            MinPar[1] = 1.0;
            MaxPar[1] = 2.1;
            NumBinsPar[1] = 24;
            LogSampling[1] = false;

            MinPar[2] = 1;
            MaxPar[2] = 5000;
            NumBinsPar[2] = 24;
            LogSampling[2] = true;
            MinPar[3] = 0.0;
            MaxPar[3] = 1.0;
            NumBinsPar[3] = 24;
            LogSampling[3] = false;

        }


    }
    else if(PotType=="SquareWell"){
        potential_id = 11;
        NumPotPars = 2;

        Kitty.SetShortRangePotential(0,0,SquareWell,cPotPars);
        //MinPar[0] = -150;
        MinPar[0] = -40;
        MaxPar[0] = -7;
        NumBinsPar[0] = 35*4;
        LogSampling[0] = false;
        //MinPar[1] = 0.4;
        MinPar[1] = 1.5;
        //MaxPar[1] = 1.6;
        MaxPar[1] = 3.5;
        NumBinsPar[1] = 44*4;
        LogSampling[1] = false;
        epsc = 5e-9;
        epsp = 1e-9;
    }
    else if(PotType=="SquareWellWall"){
        potential_id = 12;
        NumPotPars = 4;

        printf("still need to adjust pars!!!\n");
        Kitty.SetShortRangePotential(0,0,SquareWell,cPotPars);
        //MinAmp1 = -150;
        //MaxAmp1 = 0;
        //MinRng1 = 0.4;
        //MaxRng1 = 1.6;
    }

    DLM_Histo<bool> dlmTable;
    //f0,d0 and the pot pars
    dlmTable.SetUp(2+NumPotPars);
    dlmTable.SetUp(0,Num_f0, f0_goal_min, f0_goal_max);
    dlmTable.SetUp(1,Num_d0, d0_goal_min, d0_goal_max);
    for(unsigned uPar=0; uPar<NumPotPars; uPar++){
        if(NumBinsPar[uPar]==1){
            dlmTable.SetUp(2+uPar,NumBinsPar[uPar],MinPar[uPar]*0.99,MaxPar[uPar]*1.01);
        }
        else{
            dlmTable.SetUp(2+uPar,NumBinsPar[uPar],MinPar[uPar],MaxPar[uPar]);
        }
        
    }
    dlmTable.Initialize();
    dlmTable.SetBinContentAll(false);

    //Kitty.KillTheCat(); 
    Kitty.SetNotifications(CATS::nWarning);

    // before the loop
    std::ofstream* outFile = NULL;
    if(ASCII) outFile = new std::ofstream(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/output_%s_pv%i_SEED%i.txt",GetFemtoOutputFolder(), PotType.Data(), PotVar, SEED).Data());
    //std::ofstream outFile(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/output_%s_SEED%i.txt",GetFemtoOutputFolder(), PotType.Data(), SEED).Data());
    if(outFile) *outFile << "RED_MASS\tSCAT_LEN\tEFF_RNG\t\tPOT_PAR1\tPOT_PAR2\tPOT_PAR3\tPOT_PAR4\n"; 
    
    //TFile fTemp(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/froot_%s_SEED%i.root",GetFemtoOutputFolder(), PotType.Data(), SEED).Data(), "recreate");
    TFile fOutput(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/%s_%s_pv%i_SEED%i.root",GetFemtoOutputFolder(), Description.Data(), PotType.Data(), PotVar, SEED).Data(), "recreate");

    TNtuple* nt_pp_pars = new TNtuple("nt_pp_pars", "nt_pp_pars","SEED:POT_ID:RedM:f0:d0:V0:mu0:V1:mu1");
    Float_t nt_buffer[9];
    nt_buffer[0] = SEED;
    nt_buffer[1] = potential_id;


    double RedM = MinRedM;
    //a trick to have the RndPar count from zero, but we can also access the axes of our dlm_histo
    double* AxisValues = new double [6];
    double* RndPar = &AxisValues[2];
    std::vector<unsigned> evaluated_bins;

    //printf("go go %u\n", evaluated_bins.size());
    
    unsigned whileEvent = 0;
    unsigned TotBinId;
    unsigned RndBin;
    while(whileEvent<NumEvents){
        if(MinRedM!=MaxRedM) RedM = rangen.Uniform(MinRedM, MaxRedM);

        if(whileEvent%EachXevent_rnd==0 || evaluated_bins.size()==0){
            for(unsigned uPar=0; uPar<NumPotPars; uPar++){
                if(MinPar[uPar]==MaxPar[uPar]){
                    RndPar[uPar] = MinPar[uPar];
                }
                else if(LogSampling[uPar]){
                    RndPar[uPar] = rangen.Uniform(log(MinPar[uPar]), log(MaxPar[uPar]));
                    RndPar[uPar] = exp(RndPar[uPar]);                
                }
                else{
                    RndPar[uPar] = rangen.Uniform(MinPar[uPar], MaxPar[uPar]);
                }
            }
        }
        else{
            RndBin = rangen.Integer(evaluated_bins.size());
            unsigned* BinIdPerAxis = new unsigned [2+NumPotPars];
            dlmTable.GetBinCoordinates(evaluated_bins.at(RndBin), BinIdPerAxis);
            //printf("sampling from %u (%u)\n",RndBin,evaluated_bins.at(RndBin));
            for(unsigned uPar=0; uPar<NumPotPars; uPar++){
                if(MinPar[uPar]==MaxPar[uPar]){
                    RndPar[uPar] = MinPar[uPar];
                }
                else{
                    RndPar[uPar] = rangen.Uniform(dlmTable.GetBinLowEdge(2+uPar, BinIdPerAxis[2+uPar]), dlmTable.GetBinUpEdge(2+uPar, BinIdPerAxis[2+uPar]));
                }
                //printf(" %u %.3f\n", uPar, RndPar[uPar]);
                //usleep(100e3);
            }
            delete [] BinIdPerAxis;
        }
        
        /*
        //every second event is random, every second we sample from 
        //our histogram, as long as there are entries there        
        if(whileEvent%2==1 && evaluated_bins.size()>0){
            RndBin = rangen.Integer(evaluated_bins.size());
            unsigned* BinIdPerAxis = new unsigned [2+NumPotPars];
            dlmTable.GetBinCoordinates(evaluated_bins.at(RndBin), BinIdPerAxis);
            //printf("sampling from %u (%u)\n",RndBin,evaluated_bins.at(RndBin));
            for(unsigned uPar=0; uPar<NumPotPars; uPar++){
                RndPar[uPar] = rangen.Uniform(dlmTable.GetBinLowEdge(2+uPar, BinIdPerAxis[2+uPar]), dlmTable.GetBinUpEdge(2+uPar, BinIdPerAxis[2+uPar]));
                //printf(" %u %.3f\n", uPar, RndPar[uPar]);
                //usleep(500e3);
            }
            delete [] BinIdPerAxis;
        }
        else{
            for(unsigned uPar=0; uPar<NumPotPars; uPar++){
                if(LogSampling[uPar]){
                    RndPar[uPar] = rangen.Uniform(log(MinPar[uPar]), log(MaxPar[uPar]));
                    RndPar[uPar] = exp(RndPar[uPar]);                
                }
                else{
                    RndPar[uPar] = rangen.Uniform(MinPar[uPar], MaxPar[uPar]);
                }
            }
        }
        */


        /*
        if(PotType=="Gauss"){
            Amp1 = 0;
            Rng1 = 0;
            //every second event is random, every second we sample from 
            //our histogram, as long as there are entries there
            if(whileEvent%2==1){
               
            }
            if(!Amp1 && !Rng){
                Amp1 = rangen.Uniform(MinAmp1, MaxAmp1);
                Rng1 = rangen.Uniform(MinRng1, MaxRng1); 
            }

        }
        else if(PotType=="SquareWell"){
            Amp1 = rangen.Uniform(MinAmp1, MaxAmp1);
            Rng1 = rangen.Uniform(MinRng1, MaxRng1);
        }
        else if(PotType=="DoubleGauss"){
            Amp1 = rangen.Uniform(MinAmp1, MaxAmp1);
            Rng1 = rangen.Uniform(MinRng1, MaxRng1);
            Amp2 = rangen.Uniform(log(MinAmp2), log(MaxAmp2));
            Amp2 = exp(Amp2);
            Rng2 = rangen.Uniform(MinRng2, MaxRng2);    
        }
            */

        Kitty.SetRedMass(RedM);
        for(unsigned uPar=0; uPar<NumPotPars; uPar++){
            Kitty.SetShortRangePotential(0,0,uPar,RndPar[uPar]);
        }
        //Kitty.SetShortRangePotential(0,0,0,Amp1);
        //Kitty.SetShortRangePotential(0,0,1,Rng1);
        //if(PotType=="DoubleGauss"){
        //    Kitty.SetShortRangePotential(0,0,2,Amp2);
        //    Kitty.SetShortRangePotential(0,0,3,Rng2);
        //    //printf("%f %f %f %f\n",Amp1,Rng1,Amp2,Rng2);
        //}

        Kitty.SetEpsilonConv(epsc);
        Kitty.SetEpsilonProp(epsp);
        Kitty.KillTheCat();
        //printf("Done\n");

        double SCAT_LEN;
        double EFF_RNG;
        TH1F* hFIT;
        TF1* fFIT;
        if(!GetScattParameters(Kitty, SCAT_LEN, EFF_RNG, hFIT, fFIT, 2, false, false, 0)){
            continue;
        }

        //printf("%f %f %f %f\n",RndPar[0],RndPar[1],SCAT_LEN,EFF_RNG);
        //if(SCAT_LEN>17.0 && SCAT_LEN <17.6 && EFF_RNG>2.6 && EFF_RNG<3.0){
        //    printf("%f %f %f %f\n",RndPar[0],RndPar[1],SCAT_LEN,EFF_RNG);
        //    printf(" k*/tg(delta) @ k*=%.1f = %.2f\n",Kitty.GetMomentum(0),Kitty.GetMomentum(0)/tan(Kitty.GetPhaseShift(0,0,0)));
        //}

        //fTemp.cd();
        fOutput.cd();
        hFIT->SetName(TString::Format("hFIT_%i",whileEvent));
        fFIT->SetName(TString::Format("fFIT_%i",whileEvent));
        //hFIT->Write();
        //fFIT->Write();
        //printf("%p %p, %f %f\n",hFIT,fFIT,SCAT_LEN,EFF_RNG);
        delete hFIT;
        delete fFIT;

        if (std::isnan(SCAT_LEN) || std::isnan(EFF_RNG)) {
            continue;
        }
        if(fabs(EFF_RNG)>=50){
            continue;
        }

        AxisValues[0] = SCAT_LEN;
        AxisValues[1] = EFF_RNG;
        if(SCAT_LEN>=f0_goal_min && SCAT_LEN<=f0_goal_max 
            && EFF_RNG>=d0_goal_min && EFF_RNG<=d0_goal_max){
            TotBinId = dlmTable.FindTotBin(AxisValues);
            if(dlmTable.GetBinContent(TotBinId)==false){
                evaluated_bins.push_back(TotBinId);
                dlmTable.SetBinContent(TotBinId,true);
                //printf("added %.2f %.2f at %u - prs: %.2f %.2f\n",SCAT_LEN,EFF_RNG,TotBinId,RndPar[0],RndPar[1]);
                //for(unsigned uPar=0; uPar<NumPotPars+2; uPar++){
                //    printf(" av%u:%.2f ", uPar, AxisValues[uPar]);
                //}
                //printf("\n");
                //usleep(100e3);
            }
        }

        nt_buffer[2] = RedM;
        nt_buffer[3] = SCAT_LEN;
        nt_buffer[4] = EFF_RNG;

        for(unsigned uPar=0; uPar<NumPotPars; uPar++){
            nt_buffer[5+uPar] = RndPar[uPar];
        }

        fOutput.cd();
        nt_pp_pars->Fill(nt_buffer);


        if(outFile){
            *outFile << std::fixed << std::setprecision(3) << RedM << "\t\t";
            if(SCAT_LEN>=1000 || SCAT_LEN<=-100)
                *outFile << std::setprecision(2) << SCAT_LEN << "\t\t";
            else if(SCAT_LEN>=10000 || SCAT_LEN<=-1000)
                *outFile << std::setprecision(1) << SCAT_LEN << "\t\t";
            else
                *outFile << std::setprecision(3) << SCAT_LEN << "\t\t";
            if(EFF_RNG>=1000 || EFF_RNG<=-100)
                *outFile << std::setprecision(2) << EFF_RNG << "\t\t";
            else if(EFF_RNG>=10000 || EFF_RNG<=-1000)
                *outFile << std::setprecision(1) << EFF_RNG << "\t\t";
            else
                *outFile << std::setprecision(3) << EFF_RNG << "\t\t";
            *outFile << std::setprecision(3) << RndPar[0] << "\t\t";
            *outFile << std::setprecision(5) << RndPar[1] << "\t\t";
            *outFile << std::setprecision(3) << RndPar[2] << "\t\t";
            *outFile << std::setprecision(5) << RndPar[3] << "\t\t";
            *outFile << "\n";

        }

        whileEvent++;

    }

    fOutput.cd();
    nt_pp_pars->Write();
    delete nt_pp_pars;
    if(outFile) outFile->close();

    delete [] AxisValues;
    delete [] NumBinsPar;
    delete [] MinPar;
    delete [] MaxPar;
    delete [] LogSampling;
}
/*
void Generate_pp_swave_events_v1(){
    const float f0_min = 15;
    const float f0_max = 20;
    const unsigned f0_bins = 25;
    
    const float d0_min = 2;
    const float d0_max = 4;
    const unsigned d0_bins = 20;

    std::vector<std::string> pot_type;
    std::vector<std::vector<float>> par_min;
    std::vector<std::vector<float>> par_max;


    pot_type.push_back("SingleGauss");
    std::vector<float> par_min_SingleGauss;


    pot_type.push_back("DoubleGauss");
    pot_type.push_back("SquareWell");
    pot_type.push_back("SquareWellWall");



    const float d0_min = 2;
    const float d0_max = 4;
    const unsigned f0_bins = 20;

}
*/

//y= 0.0003432x2+0.05434x+3.025
//y= 0.0003432x2+0.05434x+6.025

//a function that gets the maximum number of entries in a V0 - mu0 histogram, as a function of V0.
//i.e. we get more or less a 1D function with a representation of what is V0(mu0) relation, which we can later 
//parameterize to make our life easy. This will be done for an interval of approx 30% deviation from the pp scattering pars.
//It seems the effective range has a big effect on where we pick a value on this V0(mu0), while the scattering length
//has the much more smaller effect of selecting how much up/down (d0) we deviate from the central value.
//idea: in fact one can do this for a "fixed" f0 and plot several curves showing 
//this is now offloaded to Marti

void Check_ScattPars(){

    const double kMin = 0;
    const double kMax = 60;
    const unsigned kSteps = 6;

    DLM_CommonAnaFunctions AnalysisObject; 
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

    CATS Kitty_AV18;
    Kitty_AV18.SetMomBins(kSteps, kMin, kMax);
    Kitty_AV18.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pp(Kitty_AV18,"AV18","Gauss",0,0);//McLevyNolan_Reso
    Kitty_AV18.SetQ1Q2(0);
    Kitty_AV18.KillTheCat();

    double ScatLen, EffRan;
    TH1F *hFit;
    TF1 *fitSP;
    GetScattParameters(Kitty_AV18, ScatLen, EffRan, hFit, fitSP, 2, false, false, 0);    
    printf("I got from the AV18 f0 = %.2f and d0 = %.2f\n", ScatLen, EffRan);
    printf("ERE set up (q1q2): %f\n", fitSP->GetParameter(0));

}



CATS* CatFitter = NULL;
DLM_Histo<complex<double>>* FIT_WF_U = NULL;
DLM_Histo<complex<double>>* FIT_PS = NULL;
double last_f0 = 0;
double last_d0 = 0;
const double f0_fine_lim_low = 17.3*0.9;
const double f0_fine_lim_up = 17.3*1.1;
const double d0_fine_lim_low = 2.8*0.9;
const double d0_fine_lim_up = 2.8*1.1;
    //[0] = f0
    //[1] = d0
    //[2] = norm
    //[3] = lam
    //[4] = rad
double CatFitter_pp(double* x, double* par){
    //printf("I am needed\n");
    double& kstar = *x;
    if(!CatFitter) {printf("!CatFitter\n"); return 0;}
    CatFitter->SetAnaSource(0,par[4],true);
    CatFitter->SetNotifications(CATS::nWarning);

    //f0,d0,kstar,rstar
    //void CATS::SetExternalWaveFunction(const unsigned& usCh, const unsigned& usPW, DLM_Histo<complex<double>>& histWF, DLM_Histo<complex<double>>& histPS){
    if(!FIT_WF_U){
        printf("Initial FIT_WF_U\n");

        FIT_WF_U = new DLM_Histo<complex<double>>();
        FIT_PS = new DLM_Histo<complex<double>>();

        //kstar,rstar
        FIT_WF_U->SetUp(2);
        FIT_PS->SetUp(2);

        //binning of fine and coarse should be the same in kstar and rstar. If not => there is an issue for you
        FIT_WF_U->SetUp(0,dlm_wf_u_pp_fine->GetNbins(2),dlm_wf_u_pp_fine->GetLowEdge(2),dlm_wf_u_pp_fine->GetUpEdge(2));
        FIT_PS->SetUp(0,dlm_wf_u_pp_fine->GetNbins(2),dlm_wf_u_pp_fine->GetLowEdge(2),dlm_wf_u_pp_fine->GetUpEdge(2));

        FIT_WF_U->SetUp(1,dlm_wf_u_pp_fine->GetNbins(3),dlm_wf_u_pp_fine->GetLowEdge(3),dlm_wf_u_pp_fine->GetUpEdge(3));
        FIT_PS->SetUp(1,dlm_wf_u_pp_fine->GetNbins(3),dlm_wf_u_pp_fine->GetLowEdge(3),dlm_wf_u_pp_fine->GetUpEdge(3));

        FIT_WF_U->Initialize();
        FIT_PS->Initialize();

        FIT_WF_U->SetBinContentAll(0);
        FIT_PS->SetBinContentAll(0);

        printf(" -> setup complete\n");
    }



    //update the wave function
    // dlm_wf_u_pp = f0,d0,kstar,rstar
    
    if((!last_f0&&!last_d0)||last_f0!=par[0]||last_d0!=par[1]){
        DLM_Histo<float>* dlm_wf_u_pp = NULL;
        if(par[0]>f0_fine_lim_low && par[0]<f0_fine_lim_up && par[1]>d0_fine_lim_low && par[1]<d0_fine_lim_up && dlm_wf_u_pp_fine){
            dlm_wf_u_pp = dlm_wf_u_pp_fine;
        }
        else{
            dlm_wf_u_pp = dlm_wf_u_pp_coarse;
        }

        if(!dlm_wf_u_pp) {printf("!dlm_wf_u_pp\n"); return 0;}
        unsigned WhichBin[4];
        double axis_values[4];
        WhichBin[0] = dlm_wf_u_pp->FindBin(0, par[0]);
        WhichBin[1] = dlm_wf_u_pp->FindBin(1, par[1]);
        axis_values[0] = par[0];
        axis_values[1] = par[1];
        for(unsigned uMom=0; uMom<dlm_wf_u_pp->GetNbins(2); uMom++){
            WhichBin[2] = uMom;
            axis_values[2] = dlm_wf_u_pp->GetBinCenter(2, uMom);
            for(unsigned uRad=0; uRad<dlm_wf_u_pp->GetNbins(3); uRad++){
                WhichBin[3] = uRad;
                axis_values[3] = dlm_wf_u_pp->GetBinCenter(3, uRad);
                //FIT_WF_U->SetBinContent(&WhichBin[2], dlm_wf_u_pp->GetBinContent(WhichBin));
                FIT_WF_U->SetBinContent(&WhichBin[2], dlm_wf_u_pp->Eval(axis_values));
            }
        }
        //printf("Setting up WF\n");
        CatFitter->SetExternalWaveFunction(0,0,*FIT_WF_U,*FIT_PS);
    }

    CatFitter->KillTheCat();
    double CkVal;
    CkVal = CatFitter->EvalCorrFun(kstar);
    last_f0 = par[0];
    last_d0 = par[1];
    return par[2]*(par[3]*CkVal - 1+par[3]);
}

//4 different radii
//4 different lambda pars
//attempt to make a fit with rad,lam,f0,d0 all free and see what comes out
//use s-wave only
//pot 1 Gauss, 11 is SW
void CompareToAv18(int pot_id){

    //for the cats setup
    const double kMin = 0;
    const double kMax = 400;
    const unsigned kSteps = 100;
    //for the fit
    //const double kMaxFit = 200;
    const double kMaxFit = 60;

    const double rMin = 0;
    const double rMax = 16;
    const unsigned rSteps = 320;

    TString fileName_wf_u_pp_coarse = TString::Format("%s/AI_ScatteringPars/non-AI_stuff/LookUps/wf_u_pp_Pot%i_30p.dlm",GetCernBoxDimi(),pot_id);
    TString fileName_wf_u_pp_fine = TString::Format("%s/AI_ScatteringPars/non-AI_stuff/LookUps/wf_u_pp_Pot%i_30p.dlm",GetCernBoxDimi(),pot_id);

    TFile* fInput = NULL;
    if(pot_id==11){
        fInput = new TFile("/home/dimihayl/Software/LocalFemto/Output/BaiGanio/Generate_GaussPotPars_For_Kali/SW_2025-10-30_up_to_39.root","read");
    }
    else if(pot_id==1){
        fInput = new TFile("/home/dimihayl/Software/LocalFemto/Output/BaiGanio/Generate_GaussPotPars_For_Kali/SG_2025-10-30_up_to_32.root","read");
    }

    if(!fInput){
        printf("F*ck off\n");
        return;
    }
    TNtuple* mytn = (TNtuple*)fInput->Get("nt_pp_pars");

    std::vector<float> source_size;
    source_size.push_back(1.0);
    source_size.push_back(1.5);
    source_size.push_back(2.0);
    source_size.push_back(4.0);

    std::vector<float> lambda_par;
    source_size.push_back(1.0);
    source_size.push_back(0.83);
    source_size.push_back(0.67);
    source_size.push_back(0.5);    

    DLM_CommonAnaFunctions AnalysisObject; 
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

    CATS Kitty_AV18;
    Kitty_AV18.SetMomBins(kSteps, kMin, kMax);
    //Kitty_AV18.SetNotifications(CATS::nWarning);
    //AnalysisObject.SetUpCats_pp(Kitty_AV18,"AV18_s","Gauss",0,0);
    AnalysisObject.SetUpCats_pp(Kitty_AV18,"AV18","Gauss",0,0);
    Kitty_AV18.SetEpsilonConv(2e-9);
    Kitty_AV18.SetEpsilonProp(2e-9);

    CATSparameters sPars(CATSparameters::tSource,1,true);
    CATS FatCat;
    FatCat.SetMomBins(kSteps, kMin, kMax);

    /*
    FatCat.SetAnaSource(GaussSource, sPars);
    FatCat.SetAnaSource(0, 1.0);
    FatCat.SetUseAnalyticSource(true);

    FatCat.SetMomentumDependentSource(false);
    FatCat.SetThetaDependentSource(false);
    FatCat.SetNormalizedSource(false);
    FatCat.SetAutoNormSource(true);

    FatCat.SetRedMass(0.5*Mass_p);
    FatCat.SetQ1Q2(1);
    FatCat.SetQuantumStatistics(true);
    FatCat.SetNumChannels(2);
    FatCat.SetSpin(0,0);
    FatCat.SetSpin(1,1);
    FatCat.SetNumPW(0,1);
    FatCat.SetChannelWeight(0, 0.25);
    FatCat.SetChannelWeight(1, 0.75);
    */
    AnalysisObject.SetUpCats_pp(FatCat,"AV18","Gauss",0,0);
    


    FatCat.SetNotifications(CATS::nWarning);

    //FatCat.SetEpsilonConv(2e-9);
    //FatCat.SetEpsilonProp(2e-9);

    //printf("Dont forget to make the binning finer!\n");
    
    //SetUp_Ganio_pp_pot( 96/12, 17.3*0.95, 17.3*1.05,
    //                    96/12, 2.8*0.95, 2.8*1.05,
    //                    kSteps, kMin, kMax, 
    //                    rSteps, rMin, rMax,
    //                    *mytn, pot_id);

    if(dlm_wf_u_pp_coarse) {delete dlm_wf_u_pp_coarse; dlm_wf_u_pp_coarse = NULL;}
    dlm_wf_u_pp_coarse = new DLM_Histo<float> ();
    if(!dlm_wf_u_pp_coarse->QuickLoad(fileName_wf_u_pp_coarse)){
        //printf("Loading done but bad\n");
        delete dlm_wf_u_pp_coarse;
        dlm_wf_u_pp_coarse = NULL;
        if(mytn){
            SetUp_Ganio_pp_pot( 96,17.3*0.7,17.3*1.3,
                                96,2.8*0.7,2.8*1.3,
                                kSteps, kMin, kMax, 
                                rSteps, rMin, rMax,
                                *mytn, pot_id, dlm_wf_u_pp_coarse);
            dlm_wf_u_pp_coarse->QuickWrite(TString::Format("%s/AI_ScatteringPars/non-AI_stuff/LookUps/wf_u_pp_Pot%i_30p.dlm",GetCernBoxDimi(),pot_id).Data(),true);
            
            if(dlm_wf_u_pp_fine) {delete dlm_wf_u_pp_fine; dlm_wf_u_pp_fine = NULL;}
            SetUp_Ganio_pp_pot( 96,17.3*0.9,17.3*1.1,
                                96,2.8*0.9,2.8*1.1,
                                kSteps, kMin, kMax, 
                                rSteps, rMin, rMax,
                                *mytn, pot_id, dlm_wf_u_pp_fine);   
            dlm_wf_u_pp_fine->QuickWrite(TString::Format("%s/AI_ScatteringPars/non-AI_stuff/LookUps/wf_u_pp_Pot%i_10p.dlm",GetCernBoxDimi(),pot_id).Data(),true);     
        }
        else{
            printf("I cannot load anything from anywhere\n");
            return;
        }
    }
    else{
        if(dlm_wf_u_pp_fine) {delete dlm_wf_u_pp_fine; dlm_wf_u_pp_fine = NULL;}
        dlm_wf_u_pp_fine = new DLM_Histo<float> ();
        //printf("Loading done and good\n");
        if(!dlm_wf_u_pp_fine->QuickLoad(fileName_wf_u_pp_fine)){
            printf("WARNING: dlm_wf_u_pp_fine could not be loaded!\n");
        }
        //else{
        //    printf("Fine loading done and good\n");
        //}
    }
    


    //qa for +/- 30%
    ////nt_pp_pars->Draw("d0:f0>>(96,12.11,22.49,96,1.96,3.64)","POT_ID==11","colz")
    //for a high quality fit, redo it all with +/-10%:
    //nt_pp_pars->Draw("d0:f0>>(96,15.57,19.03,96,2.52,3.08)","POT_ID==11","colz")

    

    printf("Creating the fit\n");
    CatFitter = &FatCat;
    //[0] = f0
    //[1] = d0
    //[2] = norm
    //[3] = lam
    //[4] = rad
    TF1* tf1_av18 = new TF1("tf1_av18",CatFitter_pp,kMin,kMaxFit,5);
    tf1_av18->SetParameter(0, 17.35);
    tf1_av18->SetParLimits(0, 17.3*0.9, 17.3*1.1);
    tf1_av18->SetParameter(1, 2.75);
    tf1_av18->SetParLimits(1, 2.8*0.9, 2.8*1.1);
    tf1_av18->FixParameter(2, 1);
    tf1_av18->FixParameter(3, 1);
    tf1_av18->FixParameter(4, 1.2);

    printf(" -> done\n");

    Kitty_AV18.SetAnaSource(0, tf1_av18->GetParameter(4));

    Kitty_AV18.KillTheCat();
    printf("AV18 cat killed\n");

    TFile fOutQA(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/CompareToAv18_QA.root", GetFemtoOutputFolder()),"recreate");

    //save some QA output
    //steps in f0, for several k* bins, save the WF
    for(float v_f0=17.2; v_f0<17.4; v_f0+=0.01){
        tf1_av18->SetParameter(0, v_f0);
        TGraph* gCk_qa = new TGraph();
        gCk_qa->SetName(TString::Format("gCk_qa_f%.2f",v_f0));
        for(unsigned uMom=0; uMom<CatFitter->GetNumMomBins(); uMom++){
            float v_ks = CatFitter->GetMomentum(uMom);
            if(v_ks > 80) break;
            TGraph* gwf_qa = new TGraph();
            gwf_qa->SetName(TString::Format("gwf_qa_f%.2f_k%.2f",v_f0,v_ks));
            gCk_qa->SetPoint(uMom, v_ks, tf1_av18->Eval(v_ks));
            for(float v_rad=0; v_rad < 32; v_rad += 0.05){
                gwf_qa->SetPoint(gwf_qa->GetN(), v_rad, real(CatFitter->EvalRadialWaveFunction(uMom, 0, 0, v_rad, false)));
            }
            //gCk_qa->SetPoint(gCk_qa->GetN(), v_ks, tf1_av18->Eval(v_ks));

            fOutQA.cd();
            gwf_qa->Write();
            delete gwf_qa;
        }
        fOutQA.cd();
        gCk_qa->Write();
        delete gCk_qa;
    }
    fOutQA.Close();

    
    TFile fOut(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/CompareToAv18_Pot%i.root", GetFemtoOutputFolder(),pot_id),"recreate");
    TH1F* hAv18 = new TH1F("hAv18","hAv18",kSteps, kMin, kMax);


    TGraph gAv18;
    gAv18.SetName("gAv18");
    TGraph gPotDef;
    gPotDef.SetName(TString::Format("gPotId_%i",pot_id));
    for(unsigned uKstar=0; uKstar<kSteps; uKstar++){
        double kstar, Ck_Av18, Ck_Pot;
        kstar = Kitty_AV18.GetMomentum(uKstar);
        Ck_Av18 = Kitty_AV18.GetCorrFun(uKstar);
        Ck_Pot = tf1_av18->Eval(kstar);
        //Ck_Pot = FatCat.GetCorrFun(uKstar);
        hAv18->SetBinContent(uKstar+1, Ck_Av18);
        hAv18->SetBinError(uKstar+1, 0.05);
        gAv18.SetPoint(uKstar, kstar, Ck_Av18);
        gPotDef.SetPoint(uKstar, kstar, Ck_Pot);
    }
    hAv18->Fit(tf1_av18,"S, N, R");

    TGraph gFatCatU50;
    gFatCatU50.SetName("gFatCatU50");
    TGraph gAv18U50;
    gAv18U50.SetName("gAv18U50");    

    for(double rad=0.025; rad<=8; rad+=0.05){
        gFatCatU50.SetPoint(gFatCatU50.GetN(), rad, FatCat.EvalRadialWaveFunction(FatCat.GetMomBin(50), 0, 0, rad, false).real());
        gAv18U50.SetPoint(gAv18U50.GetN(), rad, Kitty_AV18.EvalRadialWaveFunction(Kitty_AV18.GetMomBin(50), 0, 0, rad, false).real());
    }

    
    gPotDef.Write();
    gAv18.Write();
    hAv18->Write();
    tf1_av18->Write();
   
    
    gFatCatU50.Write();
    gAv18U50.Write();

    //TH2F* h_chi2_map = new TH2F("h_chi2_map","h_chi2_map",32,17.3*0.9,17.3*1.1,32,2.8*0.9,2.8*1.1);
    TH2F* h_chi2_map = new TH2F("h_chi2_map","h_chi2_map",32,17.7*0.99,17.7*1.01,32,2.75*0.99,2.75*1.01);
    for(unsigned uf0=0; uf0<h_chi2_map->GetXaxis()->GetNbins(); uf0++){
        printf("uf0 = %u\n",uf0);
        double f0_val = h_chi2_map->GetXaxis()->GetBinCenter(uf0+1);
        for(unsigned ud0=0; ud0<h_chi2_map->GetYaxis()->GetNbins(); ud0++){
            double d0_val = h_chi2_map->GetYaxis()->GetBinCenter(ud0+1);
            tf1_av18->FixParameter(0,f0_val);
            tf1_av18->FixParameter(1,d0_val);
            hAv18->Fit(tf1_av18,"Q, S, N, R");
            h_chi2_map->SetBinContent(uf0+1,ud0+1,tf1_av18->GetChisquare());
        }
    }

     h_chi2_map->Write();

    delete tf1_av18;
    delete mytn;
    delete h_chi2_map;
    delete fInput;
}

void Test_Ck_at_zero(){

    unsigned kSteps = 5;
    double kMin = 0;
    double kMax = 1;
    DLM_CommonAnaFunctions AnalysisObject; 
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

    std::vector<float> source_size;
    source_size.push_back(0.8);
    source_size.push_back(1.2);
    source_size.push_back(1.6);
    source_size.push_back(2.0);


    CATS Kitty_AV18;
    Kitty_AV18.SetMomBins(kSteps, kMin, kMax);
    //Kitty_AV18.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pp(Kitty_AV18,"AV18_s","Gauss",0,0);
    Kitty_AV18.SetQ1Q2(0);
    Kitty_AV18.SetChannelWeight(0,1);
    Kitty_AV18.SetChannelWeight(1,0);
    Kitty_AV18.SetChannelWeight(2,0);
    Kitty_AV18.SetChannelWeight(3,0);
    
    Kitty_AV18.SetEpsilonConv(1e-9);
    Kitty_AV18.SetEpsilonProp(1e-9);

    TFile fOut(TString::Format("%s/BaiGanio/Generate_GaussPotPars_For_Kali/Test_Ck_at_zero.root", GetFemtoOutputFolder()),"recreate");

    for(unsigned uSrc=0; uSrc<source_size.size(); uSrc++){
        Kitty_AV18.SetAnaSource(0,source_size.at(uSrc));
        Kitty_AV18.KillTheCat();
        TH1F* hCk = new TH1F(TString::Format("hCk_%.2f",source_size.at(uSrc)),TString::Format("hCk_%.2f",source_size.at(uSrc)),
                                kSteps, kMin, kMax);
        for(unsigned uKstar=0; uKstar<kSteps; uKstar++){
            hCk->SetBinContent(uKstar+1, Kitty_AV18.GetCorrFun(uKstar));
        }
        hCk->Write();
        delete hCk;
    }
    
}


int BAI_GANIO(int argc, char *argv[]){

//ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);  // more than default 5000
//ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3);           // or a bit looser than default 1e-4
//ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(2);             // more diagnostics

    //Test_Ck_at_zero(); return 0;
    //CompareToAv18(atoi(argv[1])); return 0;
    //Check_ScattPars(); return 0;

    //simulate enough events, so that when you run:
    //nt_pp_pars->Draw("mu0:V0>>(192,-30,-5,192,1.2,2.4)","POT_ID==1 && f0>14.3*0.998 && f0 <14.3*1.002 && d0 > 2.3*0.998 && d0<2.3*1.002","colz")
    //for each potential you get at least 150 entries. After that you make a code which takes the par values within this 0.2% interval
    //bool ASCII, TString Description, TString PotType, int PotVar, float MinRedM, float MaxRedM, unsigned NumEvents, int SEED
    Generate_GaussPotPars_For_Kali(atoi(argv[1]), argv[2], argv[3], atoi(argv[4]), Mass_p*0.5, Mass_p*0.5, atoi(argv[5]), atoi(argv[6])); return 0;
    
    //Generate_pp_swave_events_v1(argv[1], atoi(argv[2]), atoi(argv[3]));

//void SetUp_Ganio_pp_pot(const int n_f0, const double min_f0, const double max_f0,
//                        const int n_d0, const double min_d0, const double max_d0,
//                        const int n_kstar, const double min_kstar, const double max_kstar, 
//                        const int n_rstar, const double min_rstar, const double max_rstar,
//                        TNtuple& input_ntuple, const int pot_id)
    TFile fInput("/home/dimihayl/Software/LocalFemto/Output/BaiGanio/Generate_GaussPotPars_For_Kali/SW_2025-10-27_up_to_183.root","read");
    TNtuple* mytn = (TNtuple*)fInput.Get("nt_pp_pars");
    //SetUp_Ganio_pp_pot(32,17.3*0.94,17.3*1.06,32,2.8*0.94,2.8*1.06,80,0,320,320,0,16,*mytn,11);
    SetUp_Ganio_pp_pot(5,17.1,17.5,5,2.7,2.9,60,0,240,240,0,12,*mytn,2,dlm_wf_u_pp_coarse);
    //qa for +/- 30%
    ////nt_pp_pars->Draw("d0:f0>>(96,12.11,22.49,96,1.96,3.64)","POT_ID==11","colz")
    //for a high quality fit, redo it all with +/-10%:
    //nt_pp_pars->Draw("d0:f0>>(96,15.57,19.03,96,2.52,3.08)","POT_ID==11","colz")
    //TFile fInput("/home/dimihayl/Software/LocalFemto/Output/BaiGanio/Generate_GaussPotPars_For_Kali/SG_2025-10-27_up_to_165.root","read");
    //TNtuple* mytn = (TNtuple*)fInput.Get("nt_pp_pars");
    //SetUp_Ganio_pp_pot(32,17.3*0.94,17.3*1.06,32,2.8*0.94,2.8*1.06,80,0,320,320,0,16,*mytn,1);
    return 0;
}
