
#include "DimiMcPlayground.h"

#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TF1.h"

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


    TFile fOutput(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiMcPlayground/McCorrelation6D/FullRotation/fOutput_%u.root",SEED),"recreate");

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

        if(uPart_Event==0||rangen.Uniform()<0.97){
        //if(uPart_Event==0||rangen.Uniform()<1.0){
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
/*
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
*/
        }

/*
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
*/
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
        ParticleC_RE[uPart_Event].SetTheta(rangen.Uniform(0,acos(rangen.Uniform(-1.,1.))));

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

double Compute_kStar(TLorentzVector vec1, TLorentzVector vec2){
    //printf("k = %f\n",(vec1-vec2).P());
    TLorentzVector BoostVec(vec1+vec2);
    vec1.Boost(-BoostVec.BoostVector());
    vec2.Boost(-BoostVec.BoostVector());
    //printf("TotMom = %f\n",(vec1+vec2).P());
    //printf("k* = %f\n",0.5*(vec1-vec2).P());
    //printf("p1* = %f\n",(vec1).P());
    //printf("p2* = %f\n",(vec2).P());
    return 0.5*(vec1-vec2).P();
}


double Silly_pp_Ck(double& Momentum){

}

void Correlations_MomentumClasses(unsigned SEED, unsigned NumberOfEvents){

    //per event, how many classes of particles we have (with different mom distr)
    unsigned NumPart_class = 2;
    //how many event classes we have. For each one we also have a modification of the mom distr
    unsigned NumEvent_class = 3;
    const double MASS = 938;

    //!DO NOT CHANGE THIS, IT WILL NOT WORK WITH A NUMBER != 2
    const unsigned EventPartMult = 2;
    const unsigned EventPairMult = TMath::Factorial(EventPartMult-1);
    const unsigned TotParticles = NumberOfEvents*EventPartMult;
    const unsigned TotPairs = NumberOfEvents*EventPairMult;
    const unsigned MixBufferSize = EventPartMult*1;

    const unsigned MaxNumEvents_ME = 20000000;
    const unsigned TotParticles_ME = NumberOfEvents>MaxNumEvents_ME?MaxNumEvents_ME*EventPartMult:TotParticles;
    printf("TotParticles_ME = %u\n",TotParticles_ME);

    TLorentzVector* Buffer_ME_All = new TLorentzVector[TotParticles_ME];
    unsigned Num_ME_All = 0;
    TLorentzVector** Buffer_ME_Class = new TLorentzVector*[NumEvent_class];
    unsigned* Num_ME_Class = new unsigned [NumEvent_class];
    unsigned* Num_SE_Class = new unsigned [NumEvent_class];
    for(unsigned uEC=0; uEC<NumEvent_class; uEC++){
        Buffer_ME_Class[uEC] = new TLorentzVector[TotParticles_ME];
        Num_ME_Class[uEC] = 0;
        Num_SE_Class[uEC] = 0;
    }

    //with how many particles to mix each particle
    const double MixesPerParticle = 64;

    //x,y,z
    double MomMean[3][NumPart_class];
    MomMean[0][0] = 0;
    MomMean[0][1] = 0;
    MomMean[1][0] = 0;
    MomMean[1][1] = 0;
    MomMean[2][0] = 0;
    MomMean[2][1] = 0;

    //x,y,z
    double MomSigma[3][NumPart_class];
    MomSigma[0][0] = 500;
    MomSigma[0][1] = 600;
    MomSigma[1][0] = 500;
    MomSigma[1][1] = 600;
    MomSigma[2][0] = 700;
    MomSigma[2][1] = 900;

    //added to the standard momentum
    double EventMean[3][NumEvent_class];
    EventMean[0][0] = 0;
    EventMean[0][1] = 0;
    EventMean[0][2] = 0;
    EventMean[1][0] = 0;
    EventMean[1][1] = 0;
    EventMean[1][2] = 0;
    EventMean[2][0] = 0;
    EventMean[2][1] = 0;
    EventMean[2][2] = 0;

    //added to the standard momentum
    double EventSigma[3][NumEvent_class];
    EventSigma[0][0] = 0;
    EventSigma[0][1] = 100;
    EventSigma[0][2] = 200;
    EventSigma[1][0] = 0;
    EventSigma[1][1] = 100;
    EventSigma[1][2] = 200;
    EventSigma[2][0] = 0;
    EventSigma[2][1] = 100;
    EventSigma[2][2] = 200;

    TRandom3 rangen(SEED);

    TH1F* h_A = new TH1F("h_A","h_A",256,0,4096);
    TH1F* h_A_P1 = new TH1F("h_A_P1","h_A_P1",256,0,4096);
    TH1F* h_A_P2 = new TH1F("h_A_P2","h_A_P2",256,0,4096);

    TH1F* h_C = new TH1F("h_C","h_C",256,0,4096);
    TH1F* h_C_P1 = new TH1F("h_C_P1","h_C_P1",256,0,4096);
    TH1F* h_C_P2 = new TH1F("h_C_P2","h_C_P2",256,0,4096);

    TH1F* h_ME_All = new TH1F("h_ME_All","h_ME_All",256,0,4096);
    TH1F* h_ME_Class = new TH1F("h_ME_Class","h_ME_Class",256,0,4096);
    TH1F* h_ME_Class_P = new TH1F("h_ME_Class_P","h_ME_Class_P",256,0,4096);

    TH1F* h_RotPhi = new TH1F("h_RotPhi","h_RotPhi",256,0,4096);

    TH1F* hCk_CA;
    TH1F* hCk_C_ME_All;
    TH1F* hCk_C_ME_Class;
    TH1F* hCk_C_RotPhi;

    float kStar;

    for(unsigned uEvent=0; uEvent<NumberOfEvents; uEvent++){
        int EventClass = rangen.Integer(NumEvent_class);
        Num_SE_Class[EventClass] += 2;
        TLorentzVector* ParticleNonC = new TLorentzVector[EventPartMult];
        TLorentzVector* ParticleCorr = new TLorentzVector[EventPartMult];
        TLorentzVector* ParticleRotPhi = new TLorentzVector[EventPartMult];
        double RanX,RanY,RanZ;
        for(unsigned uPart=0; uPart<EventPartMult; uPart++){
            int PartClass = rangen.Integer(NumPart_class);
            RanX = rangen.Gaus(MomMean[0][PartClass]+EventMean[0][EventClass],MomSigma[0][PartClass]+EventSigma[0][EventClass]);
            RanY = rangen.Gaus(MomMean[1][PartClass]+EventMean[1][EventClass],MomSigma[1][PartClass]+EventSigma[1][EventClass]);
            RanZ = rangen.Gaus(MomMean[2][PartClass]+EventMean[2][EventClass],MomSigma[2][PartClass]+EventSigma[2][EventClass]);
//printf("MeanX = %f\n",MomMean[0][PartClass]+EventMean[0][EventClass]);
//printf("SX = %f\n",MomSigma[0][PartClass]+EventSigma[0][EventClass]);
//printf("RanX = %f\n\n",RanX);
            //Non-correlated particles
            ParticleNonC[uPart].SetXYZM(RanX,RanY,RanZ,MASS);

            //Correlated particles, where the second particle is weighted with a Fermi term
            bool Rejected = true;
            while(Rejected){
                RanX = rangen.Gaus(MomMean[0][PartClass]+EventMean[0][EventClass],MomSigma[0][PartClass]+EventSigma[0][EventClass]);
                RanY = rangen.Gaus(MomMean[1][PartClass]+EventMean[1][EventClass],MomSigma[1][PartClass]+EventSigma[1][EventClass]);
                RanZ = rangen.Gaus(MomMean[2][PartClass]+EventMean[2][EventClass],MomSigma[2][PartClass]+EventSigma[2][EventClass]);
                ParticleCorr[uPart].SetXYZM(RanX,RanY,RanZ,MASS);
                //ParticleCorr[uPart].Print();
                Rejected = false;
                if(uPart){
                    //prob for an emission in the opposite direction
                    if(rangen.Uniform()<-0.002){
                        RanX = rangen.Gaus(-ParticleCorr[0].X(),MomSigma[0][PartClass]*2.+EventSigma[0][EventClass]*2.);
                        RanY = rangen.Gaus(-ParticleCorr[0].Y(),MomSigma[1][PartClass]*2.+EventSigma[1][EventClass]*2.);
                        RanZ = rangen.Gaus(-ParticleCorr[0].Z(),MomSigma[2][PartClass]*2.+EventSigma[2][EventClass]*2.);
                        ParticleCorr[uPart].SetXYZM(RanX,RanY,RanZ,MASS);
                    }
                    else{
                        kStar = Compute_kStar(ParticleCorr[uPart],ParticleCorr[uPart-1]);
                        Rejected = (1.-0.5*exp(-pow(kStar*0.5/197.,2.0)))<rangen.Uniform();
                    }
                }
            }
        }

        if(Num_ME_All<TotParticles_ME-1){
            Buffer_ME_All[Num_ME_All] = ParticleCorr[0];
            Buffer_ME_All[Num_ME_All+1] = ParticleCorr[1];
            Num_ME_All += 2;
        }

        if(Num_ME_Class[EventClass]<TotParticles_ME-1){
            Buffer_ME_Class[EventClass][Num_ME_Class[EventClass]] = ParticleCorr[0];
            Buffer_ME_Class[EventClass][Num_ME_Class[EventClass]+1] = ParticleCorr[1];
            Num_ME_Class[EventClass] += 2;
        }

        //rotated phi
        ParticleRotPhi[0] = ParticleCorr[0];
        ParticleRotPhi[1] = ParticleCorr[1];
        ParticleRotPhi[0].SetPhi(rangen.Uniform(2.*TMath::Pi()));
        ParticleRotPhi[1].SetPhi(rangen.Uniform(2.*TMath::Pi()));
        h_RotPhi->Fill(Compute_kStar(ParticleRotPhi[0],ParticleRotPhi[1]));

        h_A->Fill(Compute_kStar(ParticleNonC[0],ParticleNonC[1]));
        h_A_P1->Fill(ParticleNonC[0].P());
        h_A_P2->Fill(ParticleNonC[1].P());

        h_C->Fill(Compute_kStar(ParticleCorr[0],ParticleCorr[1]));
        h_C_P1->Fill(ParticleCorr[0].P());
        h_C_P2->Fill(ParticleCorr[1].P());

        delete [] ParticleNonC;
        delete [] ParticleCorr;
        delete [] ParticleRotPhi;
    }

    //do the event mixing (all)
/*
    for(unsigned uPart1=0; uPart1<Num_ME_All; uPart1++){
        //+EventPartMult to avoid auto-correlation from particles from the same event
        //for(unsigned uPart2=uPart1+EventPartMult; uPart2<Num_ME_All; uPart2++){
        //    kStar = Compute_kStar(Buffer_ME_All[uPart1],Buffer_ME_All[uPart2]);
        //    h_ME_All->Fill(kStar);
        //}
        unsigned uPart2;
        for(unsigned uMix=0; uMix<MixesPerParticle; uMix++){
            uPart2 = rangen.Integer(Num_ME_All);
            kStar = Compute_kStar(Buffer_ME_All[uPart1],Buffer_ME_All[uPart2]);
            h_ME_All->Fill(kStar);
        }
    }
*/

    for(unsigned uPair=0; uPair<TotPairs; uPair++){
        unsigned uPart1 = rangen.Integer(TotParticles_ME);
        unsigned uPart2 = rangen.Integer(TotParticles_ME);
        kStar = Compute_kStar(Buffer_ME_All[uPart1],Buffer_ME_All[uPart2]);
        h_ME_All->Fill(kStar);
        //h_ME_All_P->Fill(Buffer_ME_All[uPart1].P());
        //h_ME_All_P->Fill(Buffer_ME_All[uPart2].P());
    }

    //do the event mixing (same event class)
    unsigned ParticlesToIterateOver = TotParticles_ME;
    for(unsigned uEC=0; uEC<NumEvent_class; uEC++){
        if(Num_SE_Class[uEC]<ParticlesToIterateOver) ParticlesToIterateOver = Num_SE_Class[uEC];
    }
    for(unsigned uEC=0; uEC<NumEvent_class; uEC++){
        /*
        for(unsigned uPart1=0; uPart1<ParticlesToIterateOver; uPart1++){
            unsigned uPart2;
            for(unsigned uMix=0; uMix<MixesPerParticle; uMix++){
                uPart2 = rangen.Integer(Num_ME_Class[uEC]);
                kStar = Compute_kStar(Buffer_ME_Class[uEC][uPart1],Buffer_ME_Class[uEC][uPart2]);
                h_ME_Class->Fill(kStar);
                h_ME_Class_P->Fill(Buffer_ME_Class[uEC][uPart2].P());
            }
            h_ME_Class_P->Fill(Buffer_ME_Class[uEC][uPart1].P());
        }
        */

        unsigned uPart1;
        unsigned uPart2;
        for(unsigned uPair=0; uPair<TotPairs; uPair++){
            uPart1 = rangen.Integer(Num_ME_Class[uEC]);
            uPart2 = rangen.Integer(Num_ME_Class[uEC]);
            kStar = Compute_kStar(Buffer_ME_Class[uEC][uPart1],Buffer_ME_Class[uEC][uPart2]);
            h_ME_Class->Fill(kStar);
            h_ME_Class_P->Fill(Buffer_ME_Class[uEC][uPart1].P());
            h_ME_Class_P->Fill(Buffer_ME_Class[uEC][uPart2].P());
        }

    }

    h_A->Sumw2();
    h_A_P1->Sumw2();
    h_A_P2->Sumw2();
    h_C->Sumw2();
    h_C_P1->Sumw2();
    h_C_P2->Sumw2();
    h_ME_All->Sumw2();
    h_ME_Class->Sumw2();
    h_ME_Class_P->Sumw2();
    h_RotPhi->Sumw2();

    h_A->Scale(1./h_A->Integral(),"width");
    h_A_P1->Scale(1./h_A_P1->Integral(),"width");
    h_A_P2->Scale(1./h_A_P2->Integral(),"width");
    h_C->Scale(1./h_C->Integral(),"width");
    h_C_P1->Scale(1./h_C_P1->Integral(),"width");
    h_C_P2->Scale(1./h_C_P2->Integral(),"width");
    h_ME_All->Scale(1./h_ME_All->Integral(),"width");
    h_ME_Class->Scale(1./h_ME_Class->Integral(),"width");
    h_ME_Class_P->Scale(1./h_ME_Class_P->Integral(),"width");
    h_RotPhi->Scale(1./h_RotPhi->Integral(),"width");

    hCk_CA = (TH1F*)h_C->Clone("hCk_CA");
    hCk_CA->Divide(h_A);

    hCk_C_ME_All = (TH1F*)h_C->Clone("hCk_C_ME_All");
    hCk_C_ME_All->Divide(h_ME_All);

    hCk_C_ME_Class = (TH1F*)h_C->Clone("hCk_C_ME_Class");
    hCk_C_ME_Class->Divide(h_ME_Class);

    hCk_C_RotPhi = (TH1F*)h_C->Clone("hCk_C_RotPhi");
    hCk_C_RotPhi->Divide(h_RotPhi);

    TF1* fTheoryCA = new TF1("fTheoryCA","[0]*(1-[1]*exp(-pow(x*[2]/197,2.0)))",0,4096);
    fTheoryCA->FixParameter(0,1);
    fTheoryCA->FixParameter(1,0.5);
    fTheoryCA->FixParameter(2,0.5);

    TH1F* hCkRatio_Th_A = (TH1F*)hCk_CA->Clone("hCkRatio_Th_A");
    for(unsigned uBin=1; uBin<=hCkRatio_Th_A->GetNbinsX(); uBin++){
        double kStar = hCkRatio_Th_A->GetBinCenter(uBin);
        double Ratio = fTheoryCA->Eval(kStar);
        //double RatioErr = hCk_CA->GetBinError(uBin)/fTheoryCA->Eval(kStar);
        hCkRatio_Th_A->SetBinContent(uBin,Ratio);
        hCkRatio_Th_A->SetBinError(uBin,0);
    }
    hCkRatio_Th_A->Sumw2();
    hCkRatio_Th_A->Divide(hCk_CA);

    TH1F* hCkRatio_Th_ME_Class = (TH1F*)hCk_CA->Clone("hCkRatio_Th_ME_Class");
    for(unsigned uBin=1; uBin<=hCkRatio_Th_ME_Class->GetNbinsX(); uBin++){
        double kStar = hCkRatio_Th_ME_Class->GetBinCenter(uBin);
        double Ratio = fTheoryCA->Eval(kStar);
        //double RatioErr = hCkRatio_Th_ME_Class->GetBinError(uBin)/fTheoryCA->Eval(kStar);
        hCkRatio_Th_ME_Class->SetBinContent(uBin,Ratio);
        hCkRatio_Th_ME_Class->SetBinError(uBin,0);
    }
    hCkRatio_Th_ME_Class->Sumw2();
    hCkRatio_Th_ME_Class->Divide(hCk_C_ME_Class);

    TH1F* hRatio_ME_Class_A = (TH1F*)h_ME_Class->Clone("hRatio_ME_Class_A");
    hRatio_ME_Class_A->Divide(h_A);

    TH1F* hRatio_ME_All_A = (TH1F*)h_ME_All->Clone("hRatio_ME_All_A");
    hRatio_ME_All_A->Divide(h_A);

    TH1F* hRatio_RotPhi_A = (TH1F*)h_RotPhi->Clone("hRatio_RotPhi_A");
    hRatio_RotPhi_A->Divide(h_A);

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiMcPlayground/Correlations_MomentumClasses/fOutput.root","recreate");

    h_A->Write();
    h_A_P1->Write();
    h_A_P2->Write();
    h_C->Write();
    h_C_P1->Write();
    h_C_P2->Write();
    h_ME_All->Write();
    h_ME_Class->Write();
    h_ME_Class_P->Write();
    hCk_CA->Write();
    hCk_C_ME_All->Write();
    hCk_C_ME_Class->Write();
    hCk_C_RotPhi->Write();
    fTheoryCA->Write();
    hCkRatio_Th_A->Write();
    hCkRatio_Th_ME_Class->Write();
    hRatio_ME_Class_A->Write();
    hRatio_ME_All_A->Write();
    hRatio_RotPhi_A->Write();

    delete h_A;
    delete h_A_P1;
    delete h_A_P2;

    delete h_C;
    delete h_C_P1;
    delete h_C_P2;

    delete h_ME_All;
    delete h_ME_Class;
    delete h_ME_Class_P;

    delete hCk_CA;
    delete hCk_C_ME_All;
    delete hCk_C_ME_Class;
    delete hCk_C_RotPhi;
    delete hCkRatio_Th_A;
    delete hCkRatio_Th_ME_Class;
    delete hRatio_ME_Class_A;
    delete hRatio_ME_All_A;
    delete hRatio_RotPhi_A;

    for(unsigned uEC=0; uEC<NumEvent_class; uEC++){
        delete [] Buffer_ME_Class[uEC];
    }
    delete [] Buffer_ME_Class;
    delete [] Buffer_ME_All;
    delete [] Num_ME_Class;
}



double DimiNorm_CkTh_Easy(const double& Mom, const double& Ampl, const double& Mean, const double& Stdv){
    return 1.+Ampl*TMath::Gaus(Mom,Mean,Stdv,false);
}
double DimiNorm_CkTh(double* x, double* par){
    return DimiNorm_CkTh_Easy(*x,par[0],par[1],par[2]);
}
double DimiNorm_CkTh_A(double* x, double* par){
    return (DimiNorm_CkTh(x,par)+par[3]/(1.+exp(-(*x-par[4])/par[5])))-0.5*par[3];
}
double DimiNorm_CkTh_B(double* x, double* par){
    return DimiNorm_CkTh(x,par);
}
double DimiNorm_CkTh_AB(double* x, double* par){
    return par[0]*DimiNorm_CkTh_A(x,&par[1])+(1.-par[0])*DimiNorm_CkTh_B(x,&par[7]);
}

double DimiNorm_Main(){
    const double lam_A = 0.6;

    const double mean_A = 0;
    const double stdv_A = 80;
    const double ampl_A = 3;

    const double mean_B = 100;
    const double stdv_B = 400;
    const double ampl_B = 0.6;

    const double MomentumSpread_A = 600;
    const double MomentumSpread_B = 600;
    const double MASS_A = 938;
    const double MASS_B = 938;

    const unsigned NumIter = 8000000;

    TF1* fCkTh_A = new TF1("fCkTh_A",DimiNorm_CkTh_A,0,4096,6);
    fCkTh_A->FixParameter(0,ampl_A);
    fCkTh_A->FixParameter(1,mean_A);
    fCkTh_A->FixParameter(2,stdv_A);
    fCkTh_A->FixParameter(3,0.5);
    fCkTh_A->FixParameter(4,1000);
    fCkTh_A->FixParameter(5,125);
    TF1* fCkTh_B = new TF1("fCkTh_B",DimiNorm_CkTh_B,0,4096,3);
    fCkTh_B->FixParameter(0,ampl_B);
    fCkTh_B->FixParameter(1,mean_B);
    fCkTh_B->FixParameter(2,stdv_B);
    TF1* fCkTh_AB = new TF1("fCkTh_AB",DimiNorm_CkTh_AB,0,4096,10);
    fCkTh_AB->FixParameter(0,lam_A);
    for(int iPar=0; iPar<6; iPar++) fCkTh_AB->FixParameter(1+iPar,fCkTh_A->GetParameter(iPar));
    for(int iPar=0; iPar<3; iPar++) fCkTh_AB->FixParameter(7+iPar,fCkTh_B->GetParameter(iPar));

    TH1F* hSE_A = new TH1F("hSE_A","hSE_A",1024/2,0,4096);
    TH1F* hME_A = new TH1F("hME_A","hME_A",1024/2,0,4096);
    TH1F* hSE_B = new TH1F("hSE_B","hSE_B",1024/2,0,4096);
    TH1F* hME_B = new TH1F("hME_B","hME_B",1024/2,0,4096);
    TH1F* hSE_AB = new TH1F("hSE_AB","hSE_AB",1024/2,0,4096);
    TH1F* hME_AB = new TH1F("hME_AB","hME_AB",1024/2,0,4096);

    //purity of the same event
    TH1F* hSEP_A = new TH1F("hSEP_A","hSEP_A",1024/2,0,4096);
    TH1F* hMEP_A = new TH1F("hMEP_A","hMEP_A",1024/2,0,4096);
    TH1F* hCk_PurityCorrFactor;
    TH1F* hCk_A_PurCorr;
    TH1F* hCk_B_PurCorr;


    //raw correlations
    TH1F* hCk_A;// = new TH1F("hCk_A","hCk_A",1024,0,4096);
    TH1F* hCk_B;// = new TH1F("hCk_B","hCk_B",1024,0,4096);
    TH1F* hCk_AB;// = new TH1F("hCk_AB","hCk_AB",1024,0,4096);

    //normalized at large k
    TH1F* hCk_Nlarge_A;
    TH1F* hCk_Nlarge_B;
    TH1F* hCk_Nlarge_AB;

    //normalized to the total yield
    TH1F* hCk_Nyield_A;
    TH1F* hCk_Nyield_B;
    TH1F* hCk_Nyield_AB;

    TRandom3 rangen(11);

    TLorentzVector P1_A;
    TLorentzVector P2_A;
    TLorentzVector P1_B;
    TLorentzVector P2_B;

    //TLorentzVector Pair_AA;
    //TLorentzVector Pair_BB;
    //TLorentzVector Pair_AB[4];
    double kStar_A;
    double kStar_B;
    //double kStar_AB[4];

    double momX,momY,momZ;
    double NumAcc_A=0;
    double NumAcc_B=0;

    for(unsigned uIter=0; uIter<NumIter; uIter++){
        momX = rangen.Gaus(0,MomentumSpread_A);
        momY = rangen.Gaus(0,MomentumSpread_A);
        momZ = rangen.Gaus(0,MomentumSpread_A);
        P1_A.SetXYZM(momX,momY,momZ,MASS_A);

        momX = rangen.Gaus(0,MomentumSpread_A);
        momY = rangen.Gaus(0,MomentumSpread_A);
        momZ = rangen.Gaus(0,MomentumSpread_A);
        P2_A.SetXYZM(momX,momY,momZ,MASS_A);

        momX = rangen.Gaus(0,MomentumSpread_B);
        momY = rangen.Gaus(0,MomentumSpread_B);
        momZ = rangen.Gaus(0,MomentumSpread_B);
        P1_B.SetXYZM(momX,momY,momZ,MASS_B);

        momX = rangen.Gaus(0,MomentumSpread_B);
        momY = rangen.Gaus(0,MomentumSpread_B);
        momZ = rangen.Gaus(0,MomentumSpread_B);
        P2_B.SetXYZM(momX,momY,momZ,MASS_B);

        kStar_A = Compute_kStar(P1_A,P2_A);
        kStar_B = Compute_kStar(P1_B,P2_B);

        //kStar_AB[0] = Compute_kStar(P1_A,P1_B);
        //kStar_AB[1] = Compute_kStar(P1_A,P2_B);
        //kStar_AB[2] = Compute_kStar(P2_A,P1_B);
        //kStar_AB[3] = Compute_kStar(P2_A,P2_B);

        //fill me with those guys
        hME_A->Fill(kStar_A);
        hME_B->Fill(kStar_B);
        bool TAKE_A = rangen.Uniform()<lam_A;
        if(TAKE_A){hME_AB->Fill(kStar_A);hMEP_A->Fill(kStar_A);}
        else{hME_AB->Fill(kStar_B);}

        //randomize again to have more fluctuations
        momX = rangen.Gaus(0,MomentumSpread_A);
        momY = rangen.Gaus(0,MomentumSpread_A);
        momZ = rangen.Gaus(0,MomentumSpread_A);
        P1_A.SetXYZM(momX,momY,momZ,MASS_A);

        momX = rangen.Gaus(0,MomentumSpread_A);
        momY = rangen.Gaus(0,MomentumSpread_A);
        momZ = rangen.Gaus(0,MomentumSpread_A);
        P2_A.SetXYZM(momX,momY,momZ,MASS_A);

        momX = rangen.Gaus(0,MomentumSpread_B);
        momY = rangen.Gaus(0,MomentumSpread_B);
        momZ = rangen.Gaus(0,MomentumSpread_B);
        P1_B.SetXYZM(momX,momY,momZ,MASS_B);

        momX = rangen.Gaus(0,MomentumSpread_B);
        momY = rangen.Gaus(0,MomentumSpread_B);
        momZ = rangen.Gaus(0,MomentumSpread_B);
        P2_B.SetXYZM(momX,momY,momZ,MASS_B);

        kStar_A = Compute_kStar(P1_A,P2_A);
        kStar_B = Compute_kStar(P1_B,P2_B);

        //for the se -> simply filter out
        double RAN_A = rangen.Uniform();
        double RAN_B = rangen.Uniform();
        //(ampl_A+1.) is the maximum value of the distribution
        //double ProbFor_A = DimiNorm_CkTh_Easy(kStar_A,ampl_A,mean_A,stdv_A)/(ampl_A+1.);
        //double ProbFor_A = (DimiNorm_CkTh_Easy(kStar_A,ampl_A,mean_A,stdv_A)+DimiNorm_CkTh_Easy(kStar_A,ampl_A*0.5,1800,500)-1.)/(ampl_A+1.);
        double ProbFor_A = (DimiNorm_CkTh_Easy(kStar_A,ampl_A,mean_A,stdv_A)+
                            1./(1.+exp(-(kStar_A-1000.)/150.)))
                            /(ampl_A+1.);
        double ProbFor_B = DimiNorm_CkTh_Easy(kStar_B,ampl_B,mean_B,stdv_B)/(ampl_B+1.);
        bool ACCEPT_A = ProbFor_A>RAN_A;
        bool ACCEPT_B = ProbFor_B>RAN_B;

        if(ACCEPT_A){
            hSE_A->Fill(kStar_A);
        }
        if(ACCEPT_B){
            hSE_B->Fill(kStar_B);
        }

        //if both are available, fill the histo randomly
        if(ACCEPT_A&&ACCEPT_B){
            if(TAKE_A) {hSE_AB->Fill(kStar_A);hSEP_A->Fill(kStar_A);NumAcc_A++;}
            else {hSE_AB->Fill(kStar_B);NumAcc_B++;}
        }
        //one is available, one is not
        //accept the available one, only if the fraction is biased in the other direction
        else if(ACCEPT_A&&!ACCEPT_B&&(NumAcc_A/(NumAcc_A+NumAcc_B))<lam_A){
            //hSE_AB->Fill(kStar_A);hSEP_A->Fill(kStar_A);NumAcc_A++;
        }
        else if(!ACCEPT_A&&ACCEPT_B&&(NumAcc_A/(NumAcc_A+NumAcc_B))>lam_A){
            //hSE_AB->Fill(kStar_A);NumAcc_B++;
        }
    }

    printf("The composite Ck has lam_A = %.4f (goal %.4f)\n",NumAcc_A/(NumAcc_A+NumAcc_B),lam_A);

    hSE_A->Sumw2();
    hSE_B->Sumw2();
    hSE_AB->Sumw2();
    hSEP_A->Sumw2();
    hME_A->Sumw2();
    hME_B->Sumw2();
    hME_AB->Sumw2();
    hMEP_A->Sumw2();

    hSE_A->Scale(1./hSE_A->Integral(),"width");
    hSE_B->Scale(1./hSE_B->Integral(),"width");
    hSE_AB->Scale(1./hSE_AB->Integral(),"width");
    hSEP_A->Scale(1./hSEP_A->Integral(),"width");
    hME_A->Scale(1./hME_A->Integral(),"width");
    hME_B->Scale(1./hME_B->Integral(),"width");
    hME_AB->Scale(1./hME_AB->Integral(),"width");
    hMEP_A->Scale(1./hMEP_A->Integral(),"width");

    hCk_A = (TH1F*)hSE_A->Clone("hCk_A");
    hCk_B = (TH1F*)hSE_B->Clone("hCk_B");
    hCk_AB = (TH1F*)hSE_AB->Clone("hCk_AB");

    hCk_A->Divide(hME_A);
    hCk_B->Divide(hME_B);
    hCk_AB->Divide(hME_AB);
    hSEP_A->Divide(hSE_AB);
    hMEP_A->Divide(hME_AB);

    hCk_PurityCorrFactor = (TH1F*)hSEP_A->Clone("hCk_PurityCorrFactor");
    hCk_PurityCorrFactor->Divide(hMEP_A);

    hCk_A_PurCorr = (TH1F*)hCk_AB->Clone("hCk_A_PurCorr");
    hCk_A_PurCorr->Multiply(hCk_PurityCorrFactor);

    hCk_B_PurCorr = (TH1F*)hCk_AB->Clone("hCk_B_PurCorr");
    for(unsigned uBin=0; uBin<hCk_AB->GetNbinsX(); uBin++){
        //printf("k=%.0f\n",hCk_AB->GetBinCenter(uBin+1));
        double Correction;
        Correction = 1.-lam_A*hCk_PurityCorrFactor->GetBinContent(uBin+1);
        Correction /= (1.-lam_A);
        double CorrectionError = lam_A*hCk_PurityCorrFactor->GetBinError(uBin+1);
        //printf(" err PCF = %f\n",CorrectionError);
        CorrectionError /= (1.-lam_A);
        //printf(" err CE = %f\n",CorrectionError);
        double OriginalValue = hCk_AB->GetBinContent(uBin+1);
        double OriginalError = hCk_AB->GetBinError(uBin+1);
        hCk_B_PurCorr->SetBinContent(uBin+1,OriginalValue*Correction);
        hCk_B_PurCorr->SetBinError(uBin+1,sqrt(pow(CorrectionError*OriginalError,2.)+
                                               pow(CorrectionError*OriginalValue,2.)+
                                               pow(Correction*OriginalError,2.)));
        //printf(" err OE = %f\n",OriginalError);
        //printf(" err TE = %f\n",sqrt(pow(CorrectionError*OriginalError,2.)+
        //                                       pow(CorrectionError*OriginalValue,2.)+
        //                                       pow(Correction*OriginalError,2.)));
        //hCk_B_PurCorr->SetBinError(uBin+1,OriginalError);
    }

    TFile* fOutput = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiMcPlayground/DimiNorm_Main/fOutput.root",
                               "recreate");

    fCkTh_A->Write();
    fCkTh_B->Write();
    fCkTh_AB->Write();
    hSE_A->Write();
    hME_A->Write();
    hCk_A->Write();
    hSE_B->Write();
    hME_B->Write();
    hCk_B->Write();
    hSE_AB->Write();
    hME_AB->Write();
    hCk_AB->Write();
    hSEP_A->Write();
    hMEP_A->Write();
    hCk_PurityCorrFactor->Write();
    hCk_A_PurCorr->Write();
    hCk_B_PurCorr->Write();


    delete hSE_A;
    delete hME_A;
    delete hSE_B;
    delete hME_B;
    delete hSE_AB;
    delete hME_AB;
    delete hCk_A;
    delete hCk_B;
    delete hCk_AB;
    delete hSEP_A;
    delete hMEP_A;
    delete hCk_PurityCorrFactor;
    delete hCk_A_PurCorr;
    delete hCk_B_PurCorr;
    delete fCkTh_A;
    delete fCkTh_B;
    delete fCkTh_AB;
    delete fOutput;
}


double DimiNorm2_CkTh(double *x, double* par){
    if(*x<1000) return par[0];
    else if(*x<2000) return par[1];
    else if(*x<3000) return par[2];
    else return 0;
}
double DimiNorm2_CkThSum(double *x, double* par){
    return par[0]*DimiNorm2_CkTh(x,&par[1])+(1.-par[0])*DimiNorm2_CkTh(x,&par[4]);
}

double DimiNorm2_Main(){
/*
    const double WindowLength = 1000;
    const unsigned NumWindows = 3;
    double* CkVal_A = new double [NumWindows];
    double* CkVal_B = new double [NumWindows];
    double* CkVal_AB = new double [NumWindows];
    const double lam_A = 0.5;
    TString WindowName[NumWindows+2];
    WindowName[0] = "Raw";
    WindowName[1] = "Yield";
    WindowName[2] = "W0";
    WindowName[3] = "W1";
    WindowName[4] = "W2";

    CkVal_A[0] = 3;
    CkVal_A[1] = 1;
    CkVal_A[2] = 2;

    CkVal_B[0] = 2;
    CkVal_B[1] = 2;
    CkVal_B[2] = 1;

    CkVal_AB[0] = lam_A*CkVal_A[0]+(1-lam_A)*CkVal_B[0];
    CkVal_AB[1] = lam_A*CkVal_A[1]+(1-lam_A)*CkVal_B[1];
    CkVal_AB[2] = lam_A*CkVal_A[2]+(1-lam_A)*CkVal_B[2];

    double* Weight_A = new double [NumWindows];
    double* Weight_B = new double [NumWindows];
    double Max_A=0;
    double Max_B=0;
    for(unsigned uWin=0; uWin<NumWindows; uWin++){
        if(CkVal_A[uWin]>Max_A) Max_A=CkVal_A[uWin];
        if(CkVal_B[uWin]>Max_B) Max_B=CkVal_B[uWin];
    }
    for(unsigned uWin=0; uWin<NumWindows; uWin++){
        Weight_A[uWin] = CkVal_A[uWin]/Max_A;
        Weight_B[uWin] = CkVal_B[uWin]/Max_B;
    }

    unsigned MEscale_A = 2;
    unsigned MEscale_B = 1;
    unsigned MEscale_AB = 3;

    const unsigned NumIterSE = 1000000;

    TFile* fOutput = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiMcPlayground/DimiNorm2_Main/fOutput.root",
                               "recreate");

    TF1* fCkTh_A = new TF1("fCkTh_A",DimiNorm2_CkTh,0,3000,3);
    TF1* fCkTh_B = new TF1("fCkTh_B",DimiNorm2_CkTh,0,3000,3);
    TF1* fCkTh_AB = new TF1("fCkTh_AB",DimiNorm2_CkThSum,0,3000,7);

    fCkTh_A->FixParameter(0,CkVal_A[0]);
    fCkTh_A->FixParameter(1,CkVal_A[1]);
    fCkTh_A->FixParameter(2,CkVal_A[2]);
    fCkTh_B->FixParameter(0,CkVal_B[0]);
    fCkTh_B->FixParameter(1,CkVal_B[1]);
    fCkTh_B->FixParameter(2,CkVal_B[2]);
    fCkTh_AB->FixParameter(0,lam_A);
    for(unsigned uPar=0; uPar<3; uPar++) fCkTh_AB->FixParameter(1+uPar,fCkTh_A->GetParameter(uPar));
    for(unsigned uPar=0; uPar<3; uPar++) fCkTh_AB->FixParameter(4+uPar,fCkTh_B->GetParameter(uPar));

    TH1F** hSE_A = new TH1F*[NumWindows+2];//new TH1F("hSE_A","hSE_A",512/2,0,3000);
    TH1F** hME_A = new TH1F*[NumWindows+2];//new TH1F("hME_A","hME_A",512/2,0,3000);
    TH1F** hSE_B = new TH1F*[NumWindows+2];//new TH1F("hSE_B","hSE_B",512/2,0,3000);
    TH1F** hME_B = new TH1F*[NumWindows+2];//new TH1F("hME_B","hME_B",512/2,0,3000);
    TH1F** hSE_AB = new TH1F*[NumWindows+2];//new TH1F("hSE_AB","hSE_AB",512/2,0,3000);
    TH1F** hME_AB = new TH1F*[NumWindows+2];//new TH1F("hME_AB","hME_AB",512/2,0,3000);
    for(unsigned uWin=0; uWin<NumWindows; uWin++){
        hSE_A[uWin] = new TH1F(TString::Format("hSE_%s_A",WindowName[uWin].Data()),TString::Format("hSE_%s_A",WindowName[uWin].Data()),512/2,0,3000);

    }

    //purity of the same event
    TH1F* hSEP_A = new TH1F("hSEP_A","hSEP_A",512/2,0,3000);
    TH1F* hMEP_A = new TH1F("hMEP_A","hMEP_A",512/2,0,3000);
    TH1F* hCk_PurityCorrFactor;
    TH1F** hCk_A_PurCorr = new TH1F*[NumWindows+2];
    TH1F** hCk_B_PurCorr = new TH1F*[NumWindows+2];


    //Correlations
    //normalized at
    //[0] raw
    //[1] yield
    //[2,3,4] to the window 0,1,2
    TH1F** hCk_A = new TH1F*[NumWindows+2];// = new TH1F("hCk_A","hCk_A",1024,0,4096);
    TH1F** hCk_B = new TH1F*[NumWindows+2];// = new TH1F("hCk_B","hCk_B",1024,0,4096);
    TH1F** hCk_AB = new TH1F*[NumWindows+2];// = new TH1F("hCk_AB","hCk_AB",1024,0,4096);
    //TH1F* hCk_A;
    //TH1F* hCk_B;
    //TH1F* hCk_AB;

    TRandom3 rangen(11);

    double RAND_A;
    double RAND_B;
    double RAND_AB;
    bool ACC_A;
    bool ACC_B;
    bool ACC_AB;

    for(unsigned uIter=0; uIter<NumIterSE; uIter++){
        for(unsigned uWin=0; uWin<NumWindows; uWin++){
            for(unsigned uME=0; uME<MEscale_A; uME++){
                RAND_A = rangen.Uniform(double(uWin)*1000,double(uWin)*1000+1000);
                for(unsigned uWin2=0; uWin2<NumWindows+2; uWin2++) hME_A[uWin2]->Fill(RAND_A);
            }
            for(unsigned uME=0; uME<MEscale_B; uME++){
                RAND_B = rangen.Uniform(double(uWin)*1000,double(uWin)*1000+1000);
                for(unsigned uWin2=0; uWin2<NumWindows+2; uWin2++) hME_B[uWin2]->Fill(RAND_B);
            }
            for(unsigned uME=0; uME<MEscale_AB; uME++){
                RAND_AB = rangen.Uniform(double(uWin)*1000,double(uWin)*1000+1000);
                for(unsigned uWin2=0; uWin2<NumWindows+2; uWin2++) hME_AB[uWin2]->Fill(RAND_AB);
            }

            ACC_A = rangen.Uniform()<Weight_A[uWin];
            ACC_B = rangen.Uniform()<Weight_B[uWin];
            ACC_AB = rangen.Uniform()<Weight_B[uWin];
            if(ACC_A){
                RAND_A = rangen.Uniform(double(uWin)*1000,double(uWin)*1000+1000);
                for(unsigned uWin2=0; uWin2<NumWindows+2; uWin2++) hSE_A[uWin2]->Fill(RAND_A);
            }
            if(ACC_B){
                RAND_B = rangen.Uniform(double(uWin)*1000,double(uWin)*1000+1000);
                for(unsigned uWin2=0; uWin2<NumWindows+2; uWin2++) hSE_B[uWin2]->Fill(RAND_B);
            }
            if(ACC_AB){
                RAND_AB = rangen.Uniform(double(uWin)*1000,double(uWin)*1000+1000);
                for(unsigned uWin2=0; uWin2<NumWindows+2; uWin2++) hSE_AB[uWin2]->Fill(RAND_AB);
            }
        }
    }
//hSEP_A
//hSE_A,hSE_B,hSE_AB
    //printf("The composite Ck has lam_A = %.4f (goal %.4f)\n",NumAcc_A/(NumAcc_A+NumAcc_B),lam_A);

    for(unsigned uWin=0; uWin<NumWindows+2; uWin++){
        hSE_A[uWin]->Sumw2();
        hSE_B[uWin]->Sumw2();
        hSE_AB[uWin]->Sumw2();
        hME_A[uWin]->Sumw2();
        hME_B[uWin]->Sumw2();
        hME_AB[uWin]->Sumw2();
    }

    hSE_A[1]->Scale(1./hSE_A[1]->Integral(),"width");
    hSE_A[2]->Scale(1./hSE_A[2]->Integral(1,hSE_A[2]->FindBin(1000)-1,"width"));
    hSE_A[3]->Scale(1./hSE_A[3]->Integral(hSE_A[3]->FindBin(1000)+1,hSE_A[3]->FindBin(2000)-1,"width"));
    hSE_A[4]->Scale(1./hSE_A[4]->Integral(hSE_A[4]->FindBin(2000)+1,hSE_A[4]->FindBin(3000)-1,"width"));

    hME_A[2]->Scale(1./hME_A[1]->Integral(),"width");
    hME_A[2]->Scale(1./hME_A[2]->Integral(1,hME_A[2]->FindBin(1000)-1,"width"));
    hME_A[3]->Scale(1./hME_A[3]->Integral(hME_A[3]->FindBin(1000)+1,hME_A[3]->FindBin(2000)-1,"width"));
    hME_A[4]->Scale(1./hME_A[4]->Integral(hME_A[4]->FindBin(2000)+1,hME_A[4]->FindBin(3000)-1,"width"));

    hSE_B[2]->Scale(1./hSE_B[2]->Integral(1,hSE_B[2]->FindBin(1000)-1,"width"));
    hSE_B[3]->Scale(1./hSE_B[3]->Integral(hSE_B[3]->FindBin(1000)+1,hSE_B[3]->FindBin(2000)-1,"width"));
    hSE_B[4]->Scale(1./hSE_B[4]->Integral(hSE_B[4]->FindBin(2000)+1,hSE_B[4]->FindBin(3000)-1,"width"));

    hME_B[2]->Scale(1./hME_B[2]->Integral(1,hME_B[2]->FindBin(1000)-1,"width"));
    hME_B[3]->Scale(1./hME_B[3]->Integral(hME_B[3]->FindBin(1000)+1,hME_B[3]->FindBin(2000)-1,"width"));
    hME_B[4]->Scale(1./hME_B[4]->Integral(hME_B[4]->FindBin(2000)+1,hME_B[4]->FindBin(3000)-1,"width"));

    hSE_AB[2]->Scale(1./hSE_AB[2]->Integral(1,hSE_AB[2]->FindBin(1000)-1,"width"));
    hSE_AB[3]->Scale(1./hSE_AB[3]->Integral(hSE_AB[3]->FindBin(1000)+1,hSE_AB[3]->FindBin(2000)-1,"width"));
    hSE_AB[4]->Scale(1./hSE_AB[4]->Integral(hSE_AB[4]->FindBin(2000)+1,hSE_AB[4]->FindBin(3000)-1,"width"));

    hME_AB[2]->Scale(1./hME_AB[2]->Integral(1,hME_AB[2]->FindBin(1000)-1,"width"));
    hME_AB[3]->Scale(1./hME_AB[3]->Integral(hME_AB[3]->FindBin(1000)+1,hME_AB[3]->FindBin(2000)-1,"width"));
    hME_AB[4]->Scale(1./hME_AB[4]->Integral(hME_AB[4]->FindBin(2000)+1,hME_AB[4]->FindBin(3000)-1,"width"));

    hSEP_A = (TH1F*)hSE_AB[0]->Clone("hSEP_A");
    hMEP_A = (TH1F*)hME_AB[0]->Clone("hMEP_A");

    for(unsigned uWin=0; uWin<NumWindows+2; uWin++){
        hCk_A[uWin] = (TH1F*)hSE_A[uWin]->Clone(TString::Format("hCk_%s_A",WindowName[uWin].Data(),uWin));
        hCk_B[uWin] = (TH1F*)hSE_B[uWin]->Clone(TString::Format("hCk_%s_B",WindowName[uWin].Data(),uWin));
        hCk_AB[uWin] = (TH1F*)hSE_AB[uWin]->Clone(TString::Format("hCk_%s_AB",WindowName[uWin].Data(),uWin));
    }

    hSEP_A->Sumw2();
    hMEP_A->Sumw2();
    for(unsigned uWin=0; uWin<NumWindows+2; uWin++){
        hCk_A[uWin]->Sumw2();
        hCk_B[uWin]->Sumw2();
        hCk_AB[uWin]->Sumw2();
    }

    for(unsigned uWin=0; uWin<NumWindows+2; uWin++){

    }

    hCk_A[1]->Scale(1./hCk_A->Integral(),"width");


    //hSE_A->Scale(1./hSE_A->Integral(),"width");
    //hSE_B->Scale(1./hSE_B->Integral(),"width");
    //hSE_AB->Scale(1./hSE_AB->Integral(),"width");
    //hSEP_A->Scale(1./hSEP_A->Integral(),"width");
    //hME_A->Scale(1./hME_A->Integral(),"width");
    //hME_B->Scale(1./hME_B->Integral(),"width");
    //hME_AB->Scale(1./hME_AB->Integral(),"width");
    //hMEP_A->Scale(1./hMEP_A->Integral(),"width");
    hSEP_A->Scale(hSE_AB->Integral()/hSEP_A->Integral());
    hMEP_A->Scale(hME_AB->Integral()/hMEP_A->Integral());

    for(unsigned uWin=0; uWin<NumWindows+2; uWin++){
        hCk_A[uWin]->Divide(hME_A);
        hCk_B[uWin]->Divide(hME_B);
        hCk_AB[uWin]->Divide(hME_AB);
    }
    hSEP_A->Divide(hSE_AB);
    hMEP_A->Divide(hME_AB);

    hCk_PurityCorrFactor = (TH1F*)hSEP_A->Clone("hCk_PurityCorrFactor");
    hCk_PurityCorrFactor->Divide(hMEP_A);

    for(unsigned uWin=0; uWin<NumWindows+2; uWin++){
        hCk_A_PurCorr[uWin] = (TH1F*)hCk_AB[uWin]->Clone(TString::Format("hCk_%s_A_PurCorr",WindowName.Data(),uWin));
        hCk_A_PurCorr[uWin]->Multiply(hCk_PurityCorrFactor);

        hCk_B_PurCorr[uWin] = (TH1F*)hCk_AB[uWin]->Clone(TString::Format("hCk_%s_B_PurCorr",WindowName.Data(),uWin));
        for(unsigned uBin=0; uBin<hCk_AB[uWin]->GetNbinsX(); uBin++){
            //printf("k=%.0f\n",hCk_AB->GetBinCenter(uBin+1));
            double Correction;
            Correction = 1.-lam_A*hCk_PurityCorrFactor->GetBinContent(uBin+1);
            Correction /= (1.-lam_A);
            double CorrectionError = lam_A*hCk_PurityCorrFactor->GetBinError(uBin+1);
            CorrectionError /= (1.-lam_A);
            double OriginalValue = hCk_AB[uWin]->GetBinContent(uBin+1);
            double OriginalError = hCk_AB[uWin]->GetBinError(uBin+1);
            hCk_B_PurCorr[uWin]->SetBinContent(uBin+1,OriginalValue*Correction);
            hCk_B_PurCorr[uWin]->SetBinError(uBin+1,sqrt(pow(CorrectionError*OriginalError,2.)+
                                                   pow(CorrectionError*OriginalValue,2.)+
                                                   pow(Correction*OriginalError,2.)));
        }
    }

    fCkTh_A->Write();
    fCkTh_B->Write();
    fCkTh_AB->Write();
    hSE_A->Write();
    hME_A->Write();
    for(unsigned uWin=0; uWin<NumWindows+2; uWin++) hCk_A[uWin]->Write();
    hSE_B->Write();
    hME_B->Write();
    for(unsigned uWin=0; uWin<NumWindows+2; uWin++) hCk_B[uWin]->Write();
    hSE_AB->Write();
    hME_AB->Write();
    for(unsigned uWin=0; uWin<NumWindows+2; uWin++) hCk_AB[uWin]->Write();
    hSEP_A->Write();
    hMEP_A->Write();
    hCk_PurityCorrFactor->Write();
    for(unsigned uWin=0; uWin<NumWindows+2; uWin++) hCk_A_PurCorr[uWin]->Write();
    for(unsigned uWin=0; uWin<NumWindows+2; uWin++) hCk_B_PurCorr[uWin]->Write();

    delete hSE_A;
    delete hME_A;
    delete hSE_B;
    delete hME_B;
    delete hSE_AB;
    delete hME_AB;
    for(unsigned uWin=0; uWin<NumWindows+2; uWin++){
        delete hCk_A[uWin];
        delete hCk_B[uWin];
        delete hCk_AB[uWin];
        delete hCk_A_PurCorr[uWin];
        delete hCk_B_PurCorr[uWin];
    }
    delete [] hCk_A;
    delete [] hCk_B;
    delete [] hCk_AB;
    delete hSEP_A;
    delete hMEP_A;
    delete hCk_PurityCorrFactor;
    delete [] hCk_A_PurCorr;
    delete [] hCk_B_PurCorr;
    delete fCkTh_A;
    delete fCkTh_B;
    delete fCkTh_AB;
    delete fOutput;
*/
}


int DimiMcPlayground_MAIN(int argc, char *argv[]){
    //McCorrelation6D(atoi(argv[1]),atoi(argv[2]));
    //Correlations_MomentumClasses(11,160000000*8);
    //DimiNorm_Main();
    DimiNorm2_Main();
}
