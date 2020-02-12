
#include "DimiMcPlayground.h"

#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TH2F.h"

//study what happens if you start from basics:
//independent single particle emission (A)
//independent single particle emission (B)
//a sample where I emit particles such, that each next particle depends a bit on the last (C)
//combination of particle pairs to build up C(k) for B/A and C/A
void McCorrelation6D(unsigned SEED, unsigned NumberOfEvents){
    const unsigned EventPartMult = 2;
    const unsigned EventPairMult = TMath::Factorial(EventPartMult-1);
    const unsigned TotParticles = NumberOfEvents*EventPartMult;
    const unsigned TotPairs = NumberOfEvents*EventPairMult;
    const unsigned MixBufferSize = EventPartMult*1;

    const unsigned MaxNumEvents_ME = 4000000;
    const unsigned TotParticles_ME = NumberOfEvents>MaxNumEvents_ME?MaxNumEvents_ME*EventPartMult:TotParticles;

printf("TotParticles_ME = %u\n",TotParticles_ME);
    //x,y,z
    double* MeanMom = new double [3];
    double* SigmaMom = new double [3];

    MeanMom[0] = 0;
    MeanMom[1] = 0;
    MeanMom[2] = 0;

    SigmaMom[0] = 800;
    SigmaMom[1] = 800;
    SigmaMom[2] = 1200;

    TVector3* ParticleA = new TVector3 [EventPartMult];
    TVector3* ParticleMeA = new TVector3 [TotParticles_ME];
    TVector3* ParticleMeC = new TVector3 [TotParticles_ME];

    TVector3* ParticleB = new TVector3 [EventPartMult];
    TVector3* ParticleC = new TVector3 [EventPartMult];
    TVector3* ParticleC_RE = new TVector3 [EventPartMult];
    TVector3* ParticleC_M = new TVector3 [MixBufferSize];

    TVector3* RelMom_A = new TVector3 [EventPairMult];
    TVector3* RelMom_B = new TVector3 [EventPairMult];
    TVector3* RelMom_C = new TVector3 [EventPairMult];
    TVector3* RelMom_C_RE = new TVector3 [EventPairMult];
    TVector3* RelMom_A_ME = new TVector3 [EventPairMult];
    TVector3* RelMom_C_ME = new TVector3 [EventPairMult];

    unsigned CurrentPair = 0;


    TFile* fInput = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ValeHM_ver1/CFOutput_pL_8.root","read");
    TH2F* h_pL = (TH2F*)fInput->Get("hCk_ReweightedMeV_1");
    //h_pL->Scale(1./1.8);


    TFile fOutput(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiMcPlayground/McCorrelation6D/fOutput_%u.root",SEED),"recreate");

    TH1F* hP_A = new TH1F("hP_A","hP_A",1024/2,0,4096*2);
    TH1F* hP_C = new TH1F("hP_C","hP_C",1024/2,0,4096*2);


    TH1F* hSe_A = new TH1F("hSe_A","hSe_A",1024/2,0,4096*2);
    TH1F* hSePhi_A = new TH1F("hSePhi_A","hSePhi_A",1024/2,0,6.3);
    TH1F* hSe_B = new TH1F("hSe_B","hSe_B",1024/2,0,4096*2);
    TH1F* hSe_C = new TH1F("hSe_C","hSe_C",1024/2,0,4096*2);
    TH1F* hSePhi_C = new TH1F("hSePhi_C","hSePhi_C",1024/2,0,6.3);
    TH1F* hRe_C = new TH1F("hRe_C","hRe_C",1024/2,0,4096*2);
    TH1F* hMe_A = new TH1F("hMe_A","hMe_A",1024/2,0,4096*2);
    TH1F* hMe_C = new TH1F("hMe_C","hMe_C",1024/2,0,4096*2);
    TH1F* hM_C = new TH1F("hM_C","hM_C",1024/2,0,4096);

    TH1F* hCk_AB = new TH1F("hCk_AB","hCk_AB",1024/2,0,4096*2);
    TH1F* hCk_AC = new TH1F("hCk_AC","hCk_AC",1024/2,0,4096*2);
    TH1F* hCphi_AC = new TH1F("hCphi_AC","hCphi_AC",1024/2,0,6.3);
    TH1F* hCk_ReCC = new TH1F("hCk_ReCC","hCk_ReCC",1024/2,0,4096*2);
    TH1F* hCk_ReCA = new TH1F("hCk_ReCA","hCk_ReCA",1024/2,0,4096*2);
    TH1F* hCk_MeAA = new TH1F("hCk_MeAA","hCk_MeAA",1024/2,0,4096*2);
    TH1F* hCk_MeCA = new TH1F("hCk_MeCA","hCk_MeCA",1024/2,0,4096*2);
    TH1F* hCk_MCA = new TH1F("hCk_MCA","hCk_MCA",1024/2,0,4096*2);

    TRandom3 rangen(SEED);
    unsigned NumMixedPart=0;

    for(unsigned uPart=0; uPart<TotParticles; uPart++){

        unsigned uPart_Event = uPart % EventPartMult;
//printf("uPart_Event = %u\n",uPart_Event);
        ParticleA[uPart_Event].SetX(rangen.Gaus(MeanMom[0],SigmaMom[0]));
        ParticleA[uPart_Event].SetY(rangen.Gaus(MeanMom[1],SigmaMom[1]));
        ParticleA[uPart_Event].SetZ(rangen.Gaus(MeanMom[2],SigmaMom[2]));

        ParticleB[uPart_Event].SetX(rangen.Gaus(MeanMom[0],SigmaMom[0]*1.05));
        ParticleB[uPart_Event].SetY(rangen.Gaus(MeanMom[1],SigmaMom[1]*1.05));
        ParticleB[uPart_Event].SetZ(rangen.Gaus(MeanMom[2],SigmaMom[2]*1.05));
/*
        //if(uPart_Event==0||rangen.Uniform()<0.97){
        if(uPart_Event==0||rangen.Uniform()<0.5){
            ParticleC[uPart_Event].SetX(rangen.Gaus(MeanMom[0],SigmaMom[0]));
            ParticleC[uPart_Event].SetY(rangen.Gaus(MeanMom[1],SigmaMom[1]));
            ParticleC[uPart_Event].SetZ(rangen.Gaus(MeanMom[2],SigmaMom[2]));
        }
        //correlate the emission of the particles with certain probability (3%)
        //
        else{
            ParticleC[uPart_Event].SetX(rangen.Gaus(ParticleC[uPart_Event-1].X(),SigmaMom[0]*0.15));
            ParticleC[uPart_Event].SetY(rangen.Gaus(ParticleC[uPart_Event-1].Y(),SigmaMom[1]*0.15));
            ParticleC[uPart_Event].SetZ(rangen.Gaus(ParticleC[uPart_Event-1].Z(),SigmaMom[2]*0.15));

            //induce qs
            double kStar;
            double QSprob;
            do{
                ParticleC[uPart_Event].SetX(rangen.Gaus(MeanMom[0],SigmaMom[0]));
                ParticleC[uPart_Event].SetY(rangen.Gaus(MeanMom[1],SigmaMom[1]));
                ParticleC[uPart_Event].SetZ(rangen.Gaus(MeanMom[2],SigmaMom[2]));
                kStar = (ParticleC[uPart_Event]-ParticleC[uPart_Event-1]).Mag();
                QSprob = 1-exp(-kStar*kStar/500./500.);
            }
            while(rangen.Uniform()>QSprob);

        }
*/

        //pLambda data driven correlation
        bool Rejected = true;
        while(Rejected){
            ParticleC[uPart_Event].SetX(rangen.Gaus(MeanMom[0],SigmaMom[0]));
            ParticleC[uPart_Event].SetY(rangen.Gaus(MeanMom[1],SigmaMom[1]));
            ParticleC[uPart_Event].SetZ(rangen.Gaus(MeanMom[2],SigmaMom[2]));
            double kStar;
            double dPhi;
            double REJprob;
            double CkVal;
            double CkMax = 1.8;
            double CphiMax = 0.945+0.258*exp(-pow((0-0.)/0.736,2.))+0.0109*exp(-pow((0-3.1415)/0.573,2.));
//printf("CphiMax=%f\n",CphiMax);
            const double Enhance = 1.0;
            const double Renorm = 1.0;
            if(uPart_Event==0){
                REJprob=0;
            }
            else{
                REJprob = 1;//prob to keep the particle pair
                //kStar
                kStar = (ParticleC[uPart_Event]-ParticleC[uPart_Event-1]).Mag();
                CkVal = Renorm*((h_pL->GetBinContent(h_pL->FindBin(kStar)) - 1)*Enhance + 1);
                CkMax = Renorm*((CkMax - 1)*Enhance + 1);
                //if(kStar>5900) REJprob=0.9;
                //else REJprob = 1.-(CkMax-CkVal)/CkMax;//prob to keep the particle pair
                //dPhi
                dPhi = fabs(ParticleC[uPart_Event].Phi()-ParticleC[uPart_Event-1].Phi());
                //REJprob *= TMath::Gaus(dPhi,3.1415,5.5,false);
                //tuned to reproduce experimental pAntip dPhi C
                REJprob *= (0.945+0.258*exp(-pow((dPhi-0.)/0.736,2.))+0.0109*exp(-pow((dPhi-3.1415)/0.573,2.)))/CphiMax;
                REJprob = 1.-REJprob;
            }
            Rejected = rangen.Uniform()<REJprob;
            //Rejected = false;
        }

        if(uPart<TotParticles_ME){
            ParticleMeA[uPart].SetX(ParticleA[uPart_Event].X());
            ParticleMeA[uPart].SetY(ParticleA[uPart_Event].Y());
            ParticleMeA[uPart].SetZ(ParticleA[uPart_Event].Z());

            ParticleMeC[uPart].SetX(ParticleC[uPart_Event].X());
            ParticleMeC[uPart].SetY(ParticleC[uPart_Event].Y());
            ParticleMeC[uPart].SetZ(ParticleC[uPart_Event].Z());
        }

        ParticleC_RE[uPart_Event] = ParticleC[uPart_Event];
        ParticleC_RE[uPart_Event].SetPhi(rangen.Uniform(0,2.*TMath::Pi()));
        //ParticleC_RE[uPart_Event].SetTheta(rangen.Uniform(0,acos(rangen.Uniform(-1.,1.))));

        hP_A->Fill(ParticleA[uPart_Event].Mag());
        hP_C->Fill(ParticleC[uPart_Event].Mag());

        //this is the last particle in the current event => make the pairs
        if(uPart_Event+1==EventPartMult){
            CurrentPair = 0;
            for(unsigned uP1=0; uP1<EventPartMult; uP1++){
                //same event
                for(unsigned uP2=uP1+1; uP2<EventPartMult; uP2++){
//printf(" CurrentPair = %u\n",CurrentPair);
                    RelMom_A[CurrentPair] = ParticleA[uP1]-ParticleA[uP2];
                    RelMom_B[CurrentPair] = ParticleB[uP1]-ParticleB[uP2];
                    RelMom_C[CurrentPair] = ParticleC[uP1]-ParticleC[uP2];
                    RelMom_C_RE[CurrentPair] = ParticleC_RE[uP1]-ParticleC_RE[uP2];


                    hSe_B->Fill(RelMom_B[CurrentPair].Mag());
                    hCk_AB->Fill(RelMom_B[CurrentPair].Mag());

                    hSe_C->Fill(RelMom_C[CurrentPair].Mag());
                    hSePhi_C->Fill(fabs(ParticleC[uP1].Phi()-ParticleC[uP2].Phi()));//delta phi
                    hCk_AC->Fill(RelMom_C[CurrentPair].Mag());
                    hCphi_AC->Fill(fabs(ParticleC[uP1].Phi()-ParticleC[uP2].Phi()));
                    hCk_ReCC->Fill(RelMom_C[CurrentPair].Mag());

                    hRe_C->Fill(RelMom_C_RE[CurrentPair].Mag());

                    hSe_A->Fill(RelMom_A[CurrentPair].Mag());
                    hSePhi_A->Fill(fabs(ParticleA[uP1].Phi()-ParticleA[uP2].Phi()));//delta phi
                    hCk_ReCA->Fill(RelMom_A[CurrentPair].Mag());
                    hCk_MeAA->Fill(RelMom_A[CurrentPair].Mag());
                    hCk_MeCA->Fill(RelMom_A[CurrentPair].Mag());
                    hCk_MCA->Fill(RelMom_A[CurrentPair].Mag());

                    CurrentPair++;
                }
                //mixed event
                TVector3 RelMom_C_CM;
                for(unsigned uPM=0; uPM<NumMixedPart; uPM++){
                    RelMom_C_CM = ParticleC[uP1]-ParticleC_M[uPM];
                    hM_C->Fill(RelMom_C_CM.Mag());
                }
            }

            //add stuff to the mixing buffer
            if(NumMixedPart<MixBufferSize){
                ParticleC_M[NumMixedPart] = ParticleC[uPart_Event];
                NumMixedPart++;
            }
            else{
                ParticleC_M[rangen.Integer(NumMixedPart)] = ParticleC[uPart_Event];
            }
        }
    }

    for(unsigned uEvent=0; uEvent<NumberOfEvents; uEvent++){
        TVector3* VME_A = new TVector3 [EventPartMult];
        TVector3* VME_C = new TVector3 [EventPartMult];
        for(unsigned uMult=0; uMult<EventPartMult; uMult++){
            VME_A[uMult] = ParticleMeA[rangen.Integer(TotParticles_ME)];
            VME_C[uMult] = ParticleMeC[rangen.Integer(TotParticles_ME)];
        }
        CurrentPair=0;
        for(unsigned uP1=0; uP1<EventPartMult; uP1++){
            for(unsigned uP2=uP1+1; uP2<EventPartMult; uP2++){
                RelMom_A_ME[CurrentPair] = VME_A[uP1]-VME_A[uP2];
                RelMom_C_ME[CurrentPair] = VME_C[uP1]-VME_C[uP2];
                hMe_A->Fill(RelMom_A_ME[CurrentPair].Mag());
                hMe_C->Fill(RelMom_C_ME[CurrentPair].Mag());
                CurrentPair++;
            }
        }

        delete [] VME_A;
        delete [] VME_C;
    }


    hP_A->Sumw2();
    hP_C->Sumw2();
    hSe_A->Sumw2();
    hSePhi_A->Sumw2();
    hSe_B->Sumw2();
    hSe_C->Sumw2();
    hSePhi_C->Sumw2();
    hRe_C->Sumw2();
    hMe_A->Sumw2();
    hMe_C->Sumw2();
    hM_C->Sumw2();
    hCk_AB->Sumw2();
    hCk_AC->Sumw2();
    hCphi_AC->Sumw2();
    hCk_ReCC->Sumw2();
    hCk_ReCA->Sumw2();
    hCk_MeAA->Sumw2();
    hCk_MeCA->Sumw2();
    hCk_MCA->Sumw2();

    hP_A->Scale(1./hP_A->Integral(),"");
    hP_C->Scale(1./hP_C->Integral(),"");
    hSe_A->Scale(1./hSe_A->Integral(),"");
    hSePhi_A->Scale(1./hSePhi_A->Integral(),"");
    hSe_B->Scale(1./hSe_B->Integral(),"");
    hSe_C->Scale(1./hSe_C->Integral(),"");
    hSePhi_C->Scale(1./hSePhi_C->Integral(),"");
    hRe_C->Scale(1./hRe_C->Integral(),"");
    hMe_A->Scale(1./hMe_A->Integral(),"");
    hMe_C->Scale(1./hMe_C->Integral(),"");
    hM_C->Scale(1./hM_C->Integral(),"");
    hCk_AB->Scale(1./hCk_AB->Integral(),"");
    hCk_AC->Scale(1./hCk_AC->Integral(),"");
    hCphi_AC->Scale(1./hCphi_AC->Integral(),"");
    hCk_ReCC->Scale(1./hCk_ReCC->Integral(),"");
    hCk_ReCA->Scale(1./hCk_ReCA->Integral(),"");
    hCk_MeAA->Scale(1./hCk_MeAA->Integral(),"");
    hCk_MeCA->Scale(1./hCk_MeCA->Integral(),"");
    hCk_MCA->Scale(1./hCk_MCA->Integral(),"");


    hCk_AB->Divide(hSe_A);
    hCk_AC->Divide(hSe_A);
    hCphi_AC->Divide(hSePhi_A);
    hCk_ReCC->Divide(hRe_C);
    hCk_ReCA->Divide(hRe_C);
    hCk_MeAA->Divide(hMe_A);
    hCk_MeCA->Divide(hMe_C);
    hCk_MCA->Divide(hM_C);

    hP_A->Write();
    hP_C->Write();
    hSe_A->Write();
    hSePhi_A->Write();
    hSe_B->Write();
    hSe_C->Write();
    hSePhi_C->Write();
    hRe_C->Write();
    hMe_A->Write();
    hMe_C->Write();
    hM_C->Write();
    hCk_AB->Write();
    hCk_AC->Write();
    hCphi_AC->Write();
    hCk_ReCC->Write();
    hCk_ReCA->Write();
    hCk_MeAA->Write();
    hCk_MeCA->Write();
    hCk_MCA->Write();


    delete [] MeanMom;
    delete [] SigmaMom;
    delete [] ParticleA;
    delete [] ParticleB;
    delete [] ParticleC;
    delete [] RelMom_A;
    delete [] RelMom_B;
    delete [] RelMom_C;

    delete hSe_A;
    delete hSe_B;
    delete hSe_C;
    delete hRe_C;
    delete hCk_AB;
    delete hCk_AC;
    delete hCk_ReCC;
}

int DimiMcPlayground_MAIN(int argc, char *argv[]){
    McCorrelation6D(atoi(argv[1]),atoi(argv[2]));
}
