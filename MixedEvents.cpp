#include "MixedEvents.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "DLM_CppTools.h"
#include "DLM_SubPads.h"
#include "DLM_Source.h"
#include "CommonAnaFunctions.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "DLM_Histo.h"
#include "EnvVars.h"

#include "math.h"
#include <iostream>
#include <stdio.h>

//#include "LambdaLambdaGlobalMap.h"
//#include "FemtoBoyzScripts.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TNtuple.h"

//#include "DLM_CkDecomposition.h"
//#include "DLM_CkModels.h"
//#include "DLM_Fitters.h"
//#include "DLM_Source.h"
//#include "DLM_Potentials.h"

using namespace std;

void DifferentTechniquesTest1(const TString& TranModDescr, const TString& DataSetDescr, const TString& Descr2=""){
    const double kMin = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMax = DataSetDescr=="pp"?3:DataSetDescr=="pLambda"?3:DataSetDescr=="pXim"?3:3;
    const double kMinZoom = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMaxZoom = DataSetDescr=="pp"?0.6:DataSetDescr=="pLambda"?0.6:DataSetDescr=="pXim"?0.6:0.6;
    //const double kNormMin = 0.400;
    //const double kNormMax = 0.600;
    //const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/DifferentTechniquesTest1/dEta_dPhi/";
    //const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/DifferentTechniquesTest1/ForVale/";
    const TString OutputFolder = TString::Format("%s/MixedEvents/DifferentTechniquesTest1/",GetFemtoOutputFolder());
    const unsigned NumMomBins = DataSetDescr=="pp"?100:DataSetDescr=="pLambda"?50:DataSetDescr=="pXim"?50:50;
    const TString OutFileBaseName = OutputFolder+TranModDescr+"_"+DataSetDescr+Descr2;
    TString InputFileName = DataSetDescr=="pp"?TransportFile_pp_Alice:DataSetDescr=="pLambda"?TransportFile_pL_Alice:
        DataSetDescr=="pXim"?TString("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_100KiLrz_pXim.f19"):"";
    if(DataSetDescr=="pAp"&&Descr2=="All") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_All.f19";
    if(DataSetDescr=="pAp"&&Descr2=="Core") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_Core.f19";
    if(DataSetDescr=="pAp"&&Descr2=="Corona") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_Corona.f19";
    const unsigned NumBlankHeaderLines=3;
    const unsigned MaxPairsToRead = 4e6;
    const unsigned HighMultLimit = 128;
    const unsigned MixingDepth = 8;
    const unsigned RandomDepth = 16;
    const double DimiSmearPhi = (90.)*DegToRad;
    int pdgID[2] = {2212, DataSetDescr=="pp"?2212:DataSetDescr=="pLambda"?3122:DataSetDescr=="pXim"?3312:0};
    if(DataSetDescr=="pAp"){pdgID[0]=2212; pdgID[1]=-2212;}
/*
    const unsigned NumNorm = 4;
    const double kNormMin[NumNorm] = {0.2,0.4,0.6,1};
    const double kNormMax[NumNorm] = {0.4,0.6,0.8,2};
*/
    const unsigned NumNorm = 1;
    const double kNormMin[NumNorm] = {kMin};
    const double kNormMax[NumNorm] = {kMax};

    //hNk in same event
    TH1F* hNkSe = new TH1F("hNkSe", "hNkSe", NumMomBins, kMin, kMax);
    //hNk in mixed event
    TH1F* hNkMe = new TH1F("hNkMe", "hNkMe", NumMomBins, kMin, kMax);
    //no MC intro of a signal
    TH1F* hNkSeOriginal = new TH1F("hNkSeOriginal", "hNkSeOriginal", NumMomBins, kMin, kMax);
    TH1F* hNkMeOriginal = new TH1F("hNkMeOriginal", "hNkMeOriginal", NumMomBins, kMin, kMax);

    //hNk in randomized same events
    TH1F* hNkSeRnd = new TH1F("hNkSeRnd", "hNkSeRnd", NumMomBins, kMin, kMax);
    //hNk in randomized mixed events
    TH1F* hNkMeRnd = new TH1F("hNkMeRnd", "hNkMeRnd", NumMomBins, kMin, kMax);
    //same as hNkSeRnd, but we randomize each event multiple times
    TH1F* hNkMultSeRnd = new TH1F("hNkMultSeRnd", "hNkMultSeRnd", NumMomBins, kMin, kMax);

    //hNk in DimiPhi strategy
    TH1F* hNkSeDimiPhi = new TH1F("hNkSeDimiPhi", "hNkSeDimiPhi", NumMomBins, kMin, kMax);
    //hNk in randomized mixed events
    TH1F* hNkMeDimiPhi = new TH1F("hNkMeDimiPhi", "hNkMeDimiPhi", NumMomBins, kMin, kMax);
    TH1F* hNkSeDimiPhi_Mult = new TH1F("hNkSeDimiPhi_Mult", "hNkSeDimiPhi_Mult", NumMomBins, kMin, kMax);

    TH1F* hNkMeStav = new TH1F("hNkMeStav", "hNkMeStav", NumMomBins, kMin, kMax);


    //hNkSe/hNkMe
    TH1F* hCk = new TH1F("hCk", "hCk", NumMomBins, kMin, kMax);
    TH1F* hCkOriginal = new TH1F("hCkOriginal", "hCkOriginal", NumMomBins, kMin, kMax);
    //hNkSe/hNkSeRnd
    TH1F* hCkRnd = new TH1F("hCkRnd", "hCkRnd", NumMomBins, kMin, kMax);
    //hNkSe/hNkMeRnd
    TH1F* hCkMeRnd = new TH1F("hCkMeRnd", "hCkMeRnd", NumMomBins, kMin, kMax);
    //hNkMe/hNkMeRnd
    TH1F* hCkMeMeRnd = new TH1F("hCkMeMeRnd", "hCkMeMeRnd", NumMomBins, kMin, kMax);
    //hNkSeRnd/hNkMe
    TH1F* hCkSeRndMe = new TH1F("hCkSeRndMe", "hCkSeRndMe", NumMomBins, kMin, kMax);
    //hNkSeRnd/hNkMeRnd
    TH1F* hCkSeRndMeRnd = new TH1F("hCkSeRndMeRnd", "hCkSeRndMeRnd", NumMomBins, kMin, kMax);
    //hNkSe/hNkSeDimiPhi
    TH1F* hCkDimiPhi = new TH1F("hCkDimiPhi", "hCkDimiPhi", NumMomBins, kMin, kMax);
    //hNkSe/hNkSeDimiPhi_Mult
    TH1F* hCkDimiPhi_Mult = new TH1F("hCkDimiPhi_Mult", "hCkDimiPhi_Mult", NumMomBins, kMin, kMax);
    TH1F* hCkStav = new TH1F("hCkStav", "hCkStav", NumMomBins, kMin, kMax);


    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    unsigned NumTotalPairs=0;
    //unsigned NumTotalPairsRnd=0;
    unsigned NumAllPairsOriginal=0;
    unsigned NumAllPairsRnd=0;
    unsigned NumAllPairsDimiPhi=0;
    unsigned NumAllPairsStav=0;
    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;
    unsigned uBuffer=0;
    unsigned NumSePairs;
    unsigned NumMePairs;
    unsigned NumSePairsOriginal;
    unsigned NumMePairsOriginal;
    unsigned NumSePairsRnd;
    unsigned NumSePairsDimiPhi;
    unsigned NumSePairsStav;
    //unsigned NumMePairsRnd;
    unsigned NumSePairsMultRnd;
    unsigned NumSePairsDimiPhi_Mult;

    CatsParticle KittyParticle;
    CatsParticle KittyParticleOriginal;
    CatsParticle KittyParticleRnd;
    CatsParticle* KittyParticleMultRnd = new CatsParticle[RandomDepth];

    CatsParticle KittyParticleStavinsky;
    //rotated in the phi plane (only a bit). The next particle -> rotation fixed such as to conserve E-P.
    CatsParticle KittyParticleDimiPhi;
    //same but repeated multiple times
    CatsParticle* KittyParticleDimiPhi_Mult = new CatsParticle[RandomDepth];

    CatsEvent** KittyEvent = new CatsEvent* [MixingDepth];
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        KittyEvent[uMix] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    CatsEvent** KittyEventOriginal = new CatsEvent* [MixingDepth];
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        KittyEventOriginal[uMix] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    CatsEvent** KittyEventRnd = new CatsEvent* [MixingDepth];
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        KittyEventRnd[uMix] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    //this has different structure on purpose
    CatsEvent** KittyEventMultRnd = new CatsEvent* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyEventMultRnd[uRnd] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    CatsEvent** KittyEventDimiPhi = new CatsEvent* [MixingDepth];
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        KittyEventDimiPhi[uMix] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    //this has different structure on purpose
    CatsEvent** KittyEventDimiPhi_Mult = new CatsEvent* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyEventDimiPhi_Mult[uRnd] = new CatsEvent(pdgID[0],pdgID[1]);
    }



    CatsDataBuffer* KittyBuffer = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    CatsDataBuffer* KittyBufferOriginal = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    CatsDataBuffer* KittyBufferRnd = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    CatsDataBuffer** KittyBufferMultRnd = new CatsDataBuffer* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyBufferMultRnd[uRnd] = new CatsDataBuffer(1,pdgID[0],pdgID[1]);
    }
    CatsDataBuffer* KittyBufferDimiPhi = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    CatsDataBuffer** KittyBufferDimiPhi_Mult = new CatsDataBuffer* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyBufferDimiPhi_Mult[uRnd] = new CatsDataBuffer(1,pdgID[0],pdgID[1]);
    }
    CatsLorentzVector* DimiPhi_Diff = new CatsLorentzVector[RandomDepth];

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;
TRandom3 rangen(11);
    while(!feof(InFile)){
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
//if(NumPartInEvent<1000)
//printf("NumPartInEvent=%i\n",NumPartInEvent);
        TotNumEvents++;
//printf("TotNumEvents=%u\n",TotNumEvents);
        if(NumPartInEvent>int(HighMultLimit)) RejectedHighMultEvents++;
CatsParticle FirstPart;
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            KittyParticleOriginal=KittyParticle;
            if(iPart==0) FirstPart=KittyParticle;

///induce correlation
            if(rangen.Uniform(0,1)<0.00&&iPart){
            //if(iPart){
                KittyParticle.Set(FirstPart.GetT(),FirstPart.GetX(),FirstPart.GetY(),FirstPart.GetZ(),
                                  FirstPart.GetE(),FirstPart.GetPx(),FirstPart.GetPy(),FirstPart.GetPz());
                KittyParticle.RotateMomPhi(rangen.Uniform(-20.*DegToRad,20.*DegToRad));
            }


            //KittyParticle.Set(KittyParticle.GetT(),KittyParticle.GetX(),KittyParticle.GetY(),KittyParticle.GetZ(),
            //                  KittyParticle.GetE(),KittyParticle.GetPx(),KittyParticle.GetPy(),0);
            KittyParticleRnd = KittyParticle;
            KittyParticleRnd.RotateMomPhi(rangen.Uniform(2.*TMath::Pi()));

            //KittyParticleDimi1 = KittyParticle
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                KittyParticleMultRnd[uRnd] = KittyParticle;
                KittyParticleMultRnd[uRnd].RotateMomPhi(rangen.Uniform(-DimiSmearPhi,DimiSmearPhi));

                KittyParticleDimiPhi_Mult[uRnd] = KittyParticle;
//return;
/*
                const double Smear = 0.2;
                double StretchFactor = fabs(1+rangen.Gaus(0,Smear));

                double tCrd = KittyParticleMultRnd[uRnd].GetT();
                double xCrd = KittyParticleMultRnd[uRnd].GetX();
                double yCrd = KittyParticleMultRnd[uRnd].GetY();
                double zCrd = KittyParticleMultRnd[uRnd].GetZ();
                double engy0 = KittyParticleMultRnd[uRnd].GetE();
                double xMom0 = KittyParticleMultRnd[uRnd].GetPx();
                double yMom0 = KittyParticleMultRnd[uRnd].GetPy();
                double zMom0 = KittyParticleMultRnd[uRnd].GetPz();
                double xMom=-1; xMom=xMom0*StretchFactor;
                double yMom=-1; yMom=yMom0*StretchFactor;
                double zMom=-1; zMom=zMom0*StretchFactor;
                double engy = sqrt(engy0*engy0-xMom0*xMom0-yMom0*yMom0-zMom0*zMom0+xMom*xMom+yMom*yMom+zMom*zMom);
                KittyParticleMultRnd[uRnd].Set(tCrd,xCrd,yCrd,zCrd,engy,xMom,yMom,zMom);
*/
            }
//if(NumPartInEvent>2) continue;
            if(NumTotalPairs>=MaxPairsToRead) continue;
            if(NumPartInEvent>int(HighMultLimit)) continue;

//pdgID[1] = 3122;
//if(rangen.Uniform()<0.5) KittyParticle.SetPid(3122);

//ALICE ACCEPTANCE

bool ACC_CUT = true;
if(TranModDescr=="Bhawani"&&Descr2=="with")ACC_CUT=true;
else if(TranModDescr=="Bhawani"&&Descr2=="without")ACC_CUT=false;
//printf(" ACC_CUT = %i\n",ACC_CUT);
//usleep(2000e3);

if(KittyParticle.GetP()<0.4&&ACC_CUT) continue;
if(fabs(KittyParticle.GetPseudoRap())>0.8&&ACC_CUT) continue;//!
//if(KittyParticle.GetPt()>1.5) continue;
//if(rangen.Uniform()<0.5 && KittyParticle.GetPt()<1.5) continue;

//if(KittyParticle.GetPid()==3122 && KittyParticle.GetPt()<1.5) continue;//!
//if(KittyParticle.GetPid()==3312 && KittyParticle.GetPt()<1.5) continue;//!
//if(KittyParticle.GetPid()==2212 && KittyParticle.GetP()<0.5) continue;
//if(KittyParticle.GetPid()==3122 && fabs(KittyParticle.GetPseudoRap())>1) continue;
//if(KittyParticle.GetPid()==2212 && fabs(KittyParticle.GetPseudoRap())>0.8) continue;

//printf("KittyParticle(Rnd).GetP() = %.4f (%.4f)\n",KittyParticle.GetP(),KittyParticleRnd.GetP());
//printf("KittyParticle(Rnd).GetPt() = %.4f (%.4f)\n",KittyParticle.GetPt(),KittyParticleRnd.GetPt());
//printf("KittyParticle(Rnd).GetPseudoRap() = %.4f (%.4f)\n",KittyParticle.GetPseudoRap(),KittyParticleRnd.GetPseudoRap());
//printf("\n");

//if(fabs(KittyParticle.GetPseudoRap())>0.9) printf("eta = %.2f\n",KittyParticle.GetPseudoRap());

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            KittyEvent[uBuffer]->AddParticle(KittyParticle);
            static int COUNTER=0;
//printf("wtf%i\n",COUNTER++);
//printf(" KittyParticle(Rnd).GetP() = %.4f (%.4f)\n",KittyParticle.GetP(),KittyParticleOriginal.GetP());
//printf(" KittyParticle(Rnd).GetPt() = %.4f (%.4f)\n",KittyParticle.GetPt(),KittyParticleOriginal.GetPt());
//printf(" KittyParticle(Rnd).GetPseudoRap() = %.4f (%.4f)\n",KittyParticle.GetPseudoRap(),KittyParticleOriginal.GetPseudoRap());
//printf("\n");
            KittyEventOriginal[uBuffer]->AddParticle(KittyParticleOriginal);
            KittyEventRnd[uBuffer]->AddParticle(KittyParticleRnd);
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                KittyEventMultRnd[uRnd]->AddParticle(KittyParticleMultRnd[uRnd]);
            }

            //Dimi's random phi
            //rotation are done pair-wise, such that each next rotation restores E-Mom conservation
            //the last particle is not rotated, as there are no particles left to compensate
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                if(iPart%2==0&&iPart!=NumPartInEvent){
                    KittyParticleDimiPhi_Mult[uRnd].RotateMomPhi(rangen.Uniform(-30.*DegToRad,30.*DegToRad));
                    //for(int i=0; i<100; i++)
                    //KittyParticleDimiPhi_Mult[uRnd].RotateMomPhi(rangen.Gaus(0,5.*DegToRad));

                    DimiPhi_Diff[uRnd] = KittyParticleDimiPhi_Mult[uRnd]-KittyParticle;
                }
                else if(iPart%2==1){
                    //diff = new-old
                    //new = diff+old
                    //=> new = old-diff (we want to change with -diff)
                    KittyParticleDimiPhi_Mult[uRnd] = KittyParticle-DimiPhi_Diff[uRnd];
                }
                else{
                    //nothing
                }
            }
            KittyEventDimiPhi[uBuffer]->AddParticle(KittyParticleDimiPhi_Mult[0]);
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                KittyEventDimiPhi_Mult[uRnd]->AddParticle(KittyParticleDimiPhi_Mult[uRnd]);
            }
        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)

        KittyEvent[uBuffer]->ComputeParticlePairs();
        ///KittyEventOriginal[uBuffer]->ComputeParticlePairs();
        KittyBuffer->SetEvent(uBuffer, *KittyEvent[uBuffer]);
        ///KittyBufferOriginal->SetEvent(uBuffer, *KittyEventOriginal[uBuffer]);

        KittyEventRnd[uBuffer]->ComputeParticlePairs();
        KittyBufferRnd->SetEvent(uBuffer, *KittyEventRnd[uBuffer]);

        KittyEventDimiPhi[uBuffer]->ComputeParticlePairs();
        KittyBufferDimiPhi->SetEvent(uBuffer, *KittyEventDimiPhi[uBuffer]);

        for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
            KittyEventMultRnd[uRnd]->ComputeParticlePairs();
            KittyBufferMultRnd[uRnd]->SetEvent(uBuffer, *KittyEventMultRnd[uBuffer]);

            KittyEventDimiPhi_Mult[uRnd]->ComputeParticlePairs();
            KittyBufferDimiPhi_Mult[uRnd]->SetEvent(uBuffer, *KittyEventDimiPhi_Mult[uBuffer]);
        }

        uBuffer++;

        //if the buffer is full -> empty it!
        //note that if it happens that we leave the while loop before emptying the buffer,
        //uBuffer will be != than zero! use this condition to empty the buffer when exiting the loop!
        if(uBuffer==MixingDepth){
            KittyBuffer->GoBabyGo();
            NumSePairs = KittyBuffer->GetNumPairsSameEvent();
            NumMePairs = KittyBuffer->GetNumPairsMixedEvent();
            NumTotalPairs += KittyBuffer->GetNumPairs();
            for(unsigned uSe=0; uSe<NumSePairs; uSe++){
                hNkSe->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                ///hNkSeOriginal->Fill(KittyBufferOriginal->GetSePair(uSe)->GetP()*0.5);
                hCk->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                ///hCkOriginal->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkRnd->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkMeRnd->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkDimiPhi->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkDimiPhi_Mult->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
            }
            for(unsigned uMe=0; uMe<NumMePairs; uMe++){
                hNkMe->Fill(KittyBuffer->GetMePair(uMe)->GetP()*0.5);
                hCkMeMeRnd->Fill(KittyBuffer->GetMePair(uMe)->GetP()*0.5);
            }

            ///KittyBufferOriginal->GoBabyGo();
            ///NumSePairsOriginal = KittyBufferOriginal->GetNumPairsSameEvent();
            ///NumAllPairsOriginal = KittyBufferOriginal->GetNumPairs();
            ///for(unsigned uSe=0; uSe<NumSePairsOriginal; uSe++){
            ///    hNkMeOriginal->Fill(KittyBufferOriginal->GetSePair(uSe)->GetP()*0.5);
            ///}


            KittyBufferRnd->GoBabyGo();
            NumSePairsRnd = KittyBufferRnd->GetNumPairsSameEvent();
            //NumMePairsRnd = KittyBufferRnd->GetNumPairsMixedEvent();
            NumAllPairsRnd = KittyBufferRnd->GetNumPairs();
            //double AvgNumPairsRnd = KittyBufferRnd->GetAvgNumPairs();
            for(unsigned uSe=0; uSe<NumSePairsRnd; uSe++){
                hNkSeRnd->Fill(KittyBufferRnd->GetSePair(uSe)->GetP()*0.5);
                hCkSeRndMe->Fill(KittyBufferRnd->GetSePair(uSe)->GetP()*0.5);
                hCkSeRndMeRnd->Fill(KittyBufferRnd->GetSePair(uSe)->GetP()*0.5);
            }
            //note that for the Rnd ME we can take ALL pairs, including from SE, since they should all
            //by construction carry no correlations!
            for(unsigned uMe=0; uMe<NumAllPairsRnd; uMe++){
                hNkMeRnd->Fill(KittyBufferRnd->GetPair(uMe)->GetP()*0.5);
            }

            KittyBufferDimiPhi->GoBabyGo();
            NumSePairsDimiPhi = KittyBufferDimiPhi->GetNumPairsSameEvent();
            //NumMePairsDimiPhi = KittyBufferDimiPhi->GetNumPairsMixedEvent();
//printf("NumSePairsDimiPhi=%u\n",NumSePairsDimiPhi);
            NumAllPairsDimiPhi = KittyBufferDimiPhi->GetNumPairs();
            //double AvgNumPairsDimiPhi = KittyBufferDimiPhi->GetAvgNumPairs();
            for(unsigned uSe=0; uSe<NumSePairsDimiPhi; uSe++){
                hNkSeDimiPhi->Fill(KittyBufferDimiPhi->GetSePair(uSe)->GetP()*0.5);
                //hCkSeDimiPhiMe->Fill(KittyBufferDimiPhi->GetSePair(uSe)->GetP()*0.5);
                //hCkSeDimiPhiMeDimiPhi->Fill(KittyBufferDimiPhi->GetSePair(uSe)->GetP()*0.5);
            }
            //note that for the DimiPhi ME we can take ALL pairs, including from SE, since they should all
            //by construction carry no correlations!
            for(unsigned uMe=0; uMe<NumAllPairsDimiPhi; uMe++){
                hNkMeDimiPhi->Fill(KittyBufferDimiPhi->GetPair(uMe)->GetP()*0.5);
            }

            for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
                KittyEvent[uDepth]->Reset();
                ///KittyEventOriginal[uDepth]->Reset();
                KittyEventRnd[uDepth]->Reset();
                KittyEventDimiPhi[uDepth]->Reset();
            }
            uBuffer=0;
        }

        for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
            KittyBufferMultRnd[uRnd]->GoBabyGo();
            NumSePairsMultRnd = KittyBufferMultRnd[uRnd]->GetNumPairsSameEvent();
            for(unsigned uSe=0; uSe<NumSePairsMultRnd; uSe++){
                hNkMultSeRnd->Fill(KittyBufferMultRnd[uRnd]->GetSePair(uSe)->GetP()*0.5);
            }
            KittyEventMultRnd[uRnd]->Reset();

            KittyBufferDimiPhi_Mult[uRnd]->GoBabyGo();
            NumSePairsDimiPhi_Mult = KittyBufferDimiPhi_Mult[uRnd]->GetNumPairsSameEvent();
//printf("NumSePairsDimiPhi_Mult=%u\n",NumSePairsDimiPhi_Mult);
            for(unsigned uSe=0; uSe<NumSePairsDimiPhi_Mult; uSe++){
                hNkSeDimiPhi_Mult->Fill(KittyBufferDimiPhi_Mult[uRnd]->GetSePair(uSe)->GetP()*0.5);
            }
            KittyEventDimiPhi_Mult[uRnd]->Reset();
        }

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }

    }//while(!feof(InFile))




    if(ProgressBar){
        printf("\r\033[K");
    }

    hNkSe->Sumw2();
    hNkMe->Sumw2();
    hCk->Sumw2();

    ///hNkSeOriginal->Sumw2();
    ///hNkMeOriginal->Sumw2();

    hNkSeRnd->Sumw2();
    hNkMeRnd->Sumw2();
    ///hCkOriginal->Sumw2();
    hCkRnd->Sumw2();
    hCkMeRnd->Sumw2();
    hCkMeMeRnd->Sumw2();
    hCkSeRndMe->Sumw2();
    hCkSeRndMeRnd->Sumw2();

    hNkMultSeRnd->Sumw2();

    hNkSeDimiPhi->Sumw2();
    hNkSeDimiPhi_Mult->Sumw2();
    hNkMeDimiPhi->Sumw2();
    hCkDimiPhi->Sumw2();
    hCkDimiPhi_Mult->Sumw2();

    TH1F* h_KinColl_Cons = new TH1F("h_KinColl_Cons", "h_KinColl_Cons", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        h_KinColl_Cons->SetBinContent(uBin+1, hNkMe->GetBinContent(uBin+1));
        h_KinColl_Cons->SetBinError(uBin+1, hNkMe->GetBinError(uBin+1));
    }
    h_KinColl_Cons->Divide(hNkMeRnd);
    h_KinColl_Cons->Multiply(hNkSeRnd);
    TH1F* hCk_KinColl_Cons = new TH1F("hCk_KinColl_Cons", "hCk_KinColl_Cons", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        hCk_KinColl_Cons->SetBinContent(uBin+1, hNkSe->GetBinContent(uBin+1));
        hCk_KinColl_Cons->SetBinError(uBin+1, hNkSe->GetBinError(uBin+1));
    }

    //hNkSe/hNkMultSeRnd
    TH1F* hCkMultRnd = new TH1F("hCkMultRnd", "hCkMultRnd", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        hCkMultRnd->SetBinContent(uBin+1, hNkSe->GetBinContent(uBin+1));
        hCkMultRnd->SetBinError(uBin+1, hNkSe->GetBinError(uBin+1));
    }

    //hNkSeRnd/hNkMultSeRnd
    TH1F* hCkMultRnd_SeRnd = new TH1F("hCkMultRnd_SeRnd", "hCkMultRnd_SeRnd", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        hCkMultRnd_SeRnd->SetBinContent(uBin+1, hNkSeRnd->GetBinContent(uBin+1));
        hCkMultRnd_SeRnd->SetBinError(uBin+1, hNkSeRnd->GetBinError(uBin+1));
    }

printf("hNkSe->Integral()=%e\n",hNkSe->Integral());
printf("hNkSeRnd->Integral()=%e\n",hNkSeRnd->Integral());
printf("hNkMultSeRnd->Integral()=%e\n",hNkMultSeRnd->Integral());
printf("hNkSeDimiPhi->Integral()=%e\n",hNkSeDimiPhi->Integral());
printf("hNkSeDimiPhi_Mult->Integral()=%e\n",hNkSeDimiPhi_Mult->Integral());

    hNkSe->Scale(100./hNkSe->Integral());
    hNkMe->Scale(100./hNkMe->Integral());
    hCk->Scale(100./hCk->Integral());

    ///hNkSeOriginal->Scale(100./hNkSeOriginal->Integral());
    ///hNkMeOriginal->Scale(100./hNkMeOriginal->Integral());
    ///hCkOriginal->Scale(100./hCkOriginal->Integral());
    hNkSeRnd->Scale(100./hNkSeRnd->Integral());
    hNkMeRnd->Scale(100./hNkMeRnd->Integral());
    hCkRnd->Scale(100./hCkRnd->Integral());
    hCkMeRnd->Scale(100./hCkMeRnd->Integral());
    hCkMeMeRnd->Scale(100./hCkMeMeRnd->Integral());
    hCkSeRndMe->Scale(100./hCkSeRndMe->Integral());
    hCkSeRndMeRnd->Scale(100./hCkSeRndMeRnd->Integral());
    hNkMeDimiPhi->Scale(100./hNkMeDimiPhi->Integral());
    hNkSeDimiPhi->Scale(100./hNkSeDimiPhi->Integral());
    hNkSeDimiPhi_Mult->Scale(100./hNkSeDimiPhi_Mult->Integral());
    hCkDimiPhi->Scale(100./hCkDimiPhi->Integral());
    hCkDimiPhi_Mult->Scale(100./hCkDimiPhi_Mult->Integral());

    h_KinColl_Cons->Scale(100./h_KinColl_Cons->Integral());
    hCk_KinColl_Cons->Scale(100./hCk_KinColl_Cons->Integral());

    hNkMultSeRnd->Scale(100./hNkMultSeRnd->Integral());
    hCkMultRnd->Scale(100./hCkMultRnd->Integral());
    hCkMultRnd_SeRnd->Scale(100./hCkMultRnd_SeRnd->Integral());

    hCk->Divide(hNkMe);
    ///hCkOriginal->Divide(hNkSeOriginal);
    hCkRnd->Divide(hNkSeRnd);
    hCkMeRnd->Divide(hNkMeRnd);
    hCkMeMeRnd->Divide(hNkMeRnd);
    hCkSeRndMe->Divide(hNkMe);
    hCkSeRndMeRnd->Divide(hNkMeRnd);
    hCkMultRnd->Divide(hNkMultSeRnd);
    hCkMultRnd_SeRnd->Divide(hNkMultSeRnd);

    hCkDimiPhi->Divide(hNkSeDimiPhi);
    hCkDimiPhi_Mult->Divide(hNkSeDimiPhi_Mult);

    hCk_KinColl_Cons->Divide(h_KinColl_Cons);

    //TFile* OutRoot = new TFile(TString::Format("%s.root",OutFileBaseName.Data()),"recreate");
    //OutRoot->cd();
    //hNkSe->Write();
    //hNkMe->Write();
    //hCk->Write();

    TF1* UnitLine = new TF1("UnitLine", "1", kMin, kMax);
    UnitLine->SetLineColor(kGray);
    UnitLine->SetLineWidth(5);

    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumMomBins, kMin, kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (GeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.15);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    //hAxis->GetYaxis()->SetTitle("(N per bin)x10^{-3}");
    hAxis->GetYaxis()->SetTitle("Arbitrary units");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.45);
    hAxis->GetYaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLimits(0, hNkSe->GetBinContent(hNkSe->GetMaximumBin())*1.2);
    hAxis->GetYaxis()->SetRangeUser(0, hNkSe->GetBinContent(hNkSe->GetMaximumBin())*1.2);

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", NumMomBins, kMin, kMax);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio->GetXaxis()->SetTitle("k (GeV)");
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.15);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio->GetYaxis()->SetTitle("N_{SE}/N_{ME}");
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(0.45/1.0);
    hAxisRatio->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    hAxisRatio->GetYaxis()->SetLimits(0.8, 1.2);
    hAxisRatio->GetYaxis()->SetNdivisions(504);

    TPaveText* InfoText1 = new TPaveText(0.1,0.8,0.7,0.94, "blNDC");//lbrt
    InfoText1->SetName("InfoText1");
    InfoText1->SetBorderSize(1);
    InfoText1->SetTextSize(0.065);
    InfoText1->SetFillColor(kWhite);
    InfoText1->SetTextFont(42);
    InfoText1->AddText("N.B. no multiplicity binning so far!");
    //InfoText1->AddText(TString::Format("Norm. region: [%.2f, %.2f] GeV", kNormMin, kNormMax));

    TFile* fOut = new TFile(OutFileBaseName+".root","recreate");

    for(unsigned uNorm=0; uNorm<NumNorm; uNorm++){

        TF1* NormLine = new TF1(TString::Format("NormLine%u",uNorm), "[0]", kNormMin[uNorm], kNormMax[uNorm]);
        NormLine->SetLineColor(kGray);
        NormLine->SetLineWidth(5);
        NormLine->SetParameter(0,1);

        //hAxis->GetYaxis()->SetRangeUser(0.75, 1.3);
        //hAxis->GetYaxis()->SetLimits(0.75, 1.3);


        //0 = SE/ME
        //1 = SE/RandSE
        //2 = SE/RandME
        //3 = ME/RandME
        //4 = RandSE/ME
        //5 = RandSE/RandME
        //6 = SE/SomethingThatShouldHaveOnlyCollisionKinematicsAndConservationLaws
        //7 = SE/MultSeRnd
        //8 = RandSE/MultSeRnd
        //9 = SE/DimiPhi
        //10 = SE/DimiPhi_Mult
        const unsigned NumTypes=11;

        TH1F** HistCk = new TH1F* [NumTypes];
        TH1F** HistSe = new TH1F* [NumTypes];
        TH1F** HistMe = new TH1F* [NumTypes];
        HistCk[0] = hCk;
        HistSe[0] = hNkSe;
        HistMe[0] = hNkMe;
        HistCk[1] = hCkRnd;
        HistSe[1] = hNkSe;
        HistMe[1] = hNkSeRnd;
        HistCk[2] = hCkMeRnd;
        HistSe[2] = hNkSe;
        HistMe[2] = hNkMeRnd;
        HistCk[3] = hCkMeMeRnd;
        HistSe[3] = hNkMe;
        HistMe[3] = hNkMeRnd;
        HistCk[4] = hCkSeRndMe;
        HistSe[4] = hNkSeRnd;
        HistMe[4] = hNkMe;
        HistCk[5] = hCkSeRndMeRnd;
        HistSe[5] = hNkSeRnd;
        HistMe[5] = hNkMeRnd;
        HistCk[6] = hCk_KinColl_Cons;
        HistSe[6] = hNkSe;
        HistMe[6] = h_KinColl_Cons;
        HistCk[7] = hCkMultRnd;
        HistSe[7] = hNkSe;
        HistMe[7] = hNkMultSeRnd;
        HistCk[8] = hCkMultRnd_SeRnd;
        HistSe[8] = hNkSeRnd;
        HistMe[8] = hNkMultSeRnd;
        HistCk[9] = hCkDimiPhi;
        HistSe[9] = hNkSe;
        HistMe[9] = hNkSeDimiPhi;
        HistCk[10] = hCkDimiPhi_Mult;
        HistSe[10] = hNkSe;
        HistMe[10] = hNkSeDimiPhi_Mult;
        ///HistCk[11] = hCkOriginal;
        ///HistSe[11] = hNkSe;
        ///HistMe[11] = hNkSeOriginal;

        for(unsigned uType=0; uType<NumTypes; uType++){

            HistCk[uType]->Fit(NormLine, "S, N, R, M");
            HistSe[uType]->Scale(1./NormLine->GetParameter(0));
            HistCk[uType]->Scale(1./NormLine->GetParameter(0));

            HistSe[uType]->SetLineWidth(6);
            HistSe[uType]->SetLineColor(kBlue+1);

            HistMe[uType]->SetLineWidth(6);
            HistMe[uType]->SetLineColor(kRed+1);

            HistCk[uType]->SetLineWidth(6);
            HistCk[uType]->SetLineColor(kBlack);

            TH1F* hShade = new TH1F(TString::Format("hShade%u_%u",uNorm,uType), TString::Format("hShade%u_%u",uNorm,uType), 1, kNormMin[uNorm], kNormMax[uNorm]);
            double hShadeMaxVal = HistSe[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2;
            if(HistCk[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2>hShadeMaxVal) hShadeMaxVal = HistCk[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2;
            hShade->SetBinContent(1, hShadeMaxVal);
            hShade->SetFillStyle(3344);//4050
            hShade->SetFillColor(46);
            hShade->SetLineColor(46);
            hShade->SetLineWidth(0);

            TLegend* myLegend = new TLegend(0.7,0.5,0.96,0.94);//lbrt
            myLegend->SetName(TString::Format("myLegend%u_%u",uNorm,uType));
            myLegend->SetTextSize(0.07);
            if(uType==0){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS mixed event");
            }
            else if(uType==1){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized SE");
            }
            else if(uType==2){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized ME");
            }
            else if(uType==3){
                myLegend->AddEntry(HistSe[uType], "EPOS mixed event");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized ME");
            }
            else if(uType==4){
                myLegend->AddEntry(HistSe[uType], "EPOS randomized SE");
                myLegend->AddEntry(HistMe[uType], "EPOS ME");
            }
            else if(uType==5){
                myLegend->AddEntry(HistSe[uType], "EPOS randomized SE");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized ME");
            }
            else if(uType==6){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Dimi's fantasy");
            }
            else if(uType==7){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Andi's fantasy");
            }
            else if(uType==8){
                myLegend->AddEntry(HistSe[uType], "EPOS randomized SE");
                myLegend->AddEntry(HistMe[uType], "Andi's fantasy");
            }
            else if(uType==9){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Dimi random phi");
            }
            else if(uType==10){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Dimi random phi (mult)");
            }
            else if(uType==11){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS original event");
            }

            myLegend->AddEntry(HistCk[uType], "Ratio");
            myLegend->AddEntry(hShade, "Norm. region");

            DLM_SubPads DrawBoard(1920,1080);
            DrawBoard.AddSubPadTL(0,1,0,0.5);//lrtb
            DrawBoard.AddSubPadTL(0,1,0.5,1);//lrtb

            DrawBoard.SetMargin(0, 0.1, 0.04, 0.0, 0.05);//lrbt
            DrawBoard.SetMargin(1, 0.1, 0.04, 0.15, 0);//lrbt

        ////////////
            DrawBoard.cd(0);

            DrawBoard.SetLabelSize(0, hAxis->GetXaxis(), 18);
            DrawBoard.SetLabelSize(0, hAxis->GetYaxis(), 18);
            DrawBoard.SetTitleSize(0, hAxis->GetXaxis(), 18);
            DrawBoard.SetTitleSize(0, hAxis->GetYaxis(), 18);

            hAxis->GetYaxis()->SetRangeUser(0, HistSe[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2);
            hAxis->GetXaxis()->SetRangeUser(kMin,kMax);
            hAxisRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
            hAxisRatio->GetXaxis()->SetRangeUser(kMin,kMax);

            hAxis->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            HistSe[uType]->Draw("same");
            HistMe[uType]->Draw("same");
            myLegend->Draw("same");
            InfoText1->Draw("same");

            DrawBoard.cd(1);

            DrawBoard.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoard.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
            DrawBoard.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoard.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

            hAxisRatio->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            UnitLine->Draw("same");
            HistCk[uType]->Draw("same");

            Width_t width_hShade=hShade->GetLineWidth();
            Width_t width_hNkSe=HistSe[uType]->GetLineWidth();
            Width_t width_hNkMe=HistMe[uType]->GetLineWidth();
            Width_t width_UnitLine=UnitLine->GetLineWidth();
            Width_t width_hCk=HistCk[uType]->GetLineWidth();

            DrawBoard.GetCanvas()->SaveAs(TString::Format("%s_DB_N%.0f-%.0f_%u.png",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade/2.5);
            HistSe[uType]->SetLineWidth(width_hNkSe/2.5);
            HistMe[uType]->SetLineWidth(width_hNkMe/2.5);
            UnitLine->SetLineWidth(width_UnitLine/2.5);
            HistCk[uType]->SetLineWidth(width_hCk/2.5);
            DrawBoard.GetCanvas()->SaveAs(TString::Format("%s_DB_N%.0f-%.0f_%u.pdf",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade);
            HistSe[uType]->SetLineWidth(width_hNkSe);
            HistMe[uType]->SetLineWidth(width_hNkMe);
            UnitLine->SetLineWidth(width_UnitLine);
            HistCk[uType]->SetLineWidth(width_hCk);

        /////////////////////////////////////////////////////////////////////
            DLM_SubPads DrawBoardZoom(1920,1080);
            DrawBoardZoom.AddSubPadTL(0,1,0,0.5);//lrtb
            DrawBoardZoom.AddSubPadTL(0,1,0.5,1);//lrtb

            DrawBoardZoom.SetMargin(0, 0.1, 0.04, 0.0, 0.05);//lrbt
            DrawBoardZoom.SetMargin(1, 0.1, 0.04, 0.15, 0);//lrbt

        ////////////
            DrawBoardZoom.cd(0);

            DrawBoardZoom.SetLabelSize(0, hAxis->GetXaxis(), 18);
            DrawBoardZoom.SetLabelSize(0, hAxis->GetYaxis(), 18);
            DrawBoardZoom.SetTitleSize(0, hAxis->GetXaxis(), 18);
            DrawBoardZoom.SetTitleSize(0, hAxis->GetYaxis(), 18);

            unsigned BinLow = HistSe[uType]->FindBin(kMinZoom);
            unsigned BinUp = HistSe[uType]->FindBin(kMaxZoom);
            double MaxYval=0;
            for(unsigned uBin=BinLow; uBin<=BinUp; uBin++){
                if(MaxYval<HistSe[uType]->GetBinContent(uBin)) MaxYval=HistSe[uType]->GetBinContent(uBin);
            }
            hAxis->GetYaxis()->SetRangeUser(0, MaxYval*1.2);
            hAxis->GetYaxis()->SetLimits(0, MaxYval*1.2);
            hAxis->GetXaxis()->SetRangeUser(kMinZoom, kMaxZoom);
            //hAxis->GetXaxis()->SetLimits(kMinZoom, kMaxZoom);

            hAxisRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
            hAxisRatio->GetYaxis()->SetLimits(0.8, 1.2);
            hAxisRatio->GetXaxis()->SetRangeUser(kMinZoom, kMaxZoom);
            //hAxisRatio->GetXaxis()->SetLimits(kMinZoom, kMaxZoom);

            hAxis->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            HistSe[uType]->Draw("same");
            HistMe[uType]->Draw("same");
            myLegend->Draw("same");
            InfoText1->Draw("same");

            DrawBoardZoom.cd(1);

            DrawBoardZoom.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoardZoom.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
            DrawBoardZoom.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoardZoom.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

            hAxisRatio->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            UnitLine->Draw("same");
            HistCk[uType]->Draw("same");

            DrawBoardZoom.GetCanvas()->SaveAs(TString::Format("%s_DBZoom_N%.0f-%.0f_%u.png",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade/2.5);
            HistSe[uType]->SetLineWidth(width_hNkSe/2.5);
            HistMe[uType]->SetLineWidth(width_hNkMe/2.5);
            UnitLine->SetLineWidth(width_UnitLine/2.5);
            HistCk[uType]->SetLineWidth(width_hCk/2.5);
            DrawBoardZoom.GetCanvas()->SaveAs(TString::Format("%s_DBZoom_N%.0f-%.0f_%u.pdf",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade);
            HistSe[uType]->SetLineWidth(width_hNkSe);
            HistMe[uType]->SetLineWidth(width_hNkMe);
            UnitLine->SetLineWidth(width_UnitLine);
            HistCk[uType]->SetLineWidth(width_hCk);

            HistSe[uType]->Write();
            HistMe[uType]->Write();
            HistCk[uType]->Write();

            delete hShade;
            delete myLegend;
        }

        delete [] HistCk;
        delete [] HistSe;
        delete [] HistMe;

        delete NormLine;

    }

    delete hNkSe;
    delete hNkSeOriginal;
    delete hNkSeRnd;
    delete hNkMeRnd;
    delete hNkMe;
    delete hNkMeOriginal;
    delete hNkMeDimiPhi;
    delete hNkSeDimiPhi;
    delete hNkSeDimiPhi_Mult;
    delete hCk;
    delete hCkOriginal;
    delete hCkRnd;
    delete hCkMeRnd;
    delete hCkMeMeRnd;
    delete hCkSeRndMe;
    delete hCkSeRndMeRnd;
    delete hCkMultRnd;
    delete hCkMultRnd_SeRnd;
    delete hCkDimiPhi;
    delete hCkDimiPhi_Mult;

    delete h_KinColl_Cons;
    delete hCk_KinColl_Cons;

    delete UnitLine;
    delete hAxis;
    delete hAxisRatio;
    delete InfoText1;
    //delete OutRoot;

    delete [] cdummy;
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        delete KittyEvent[uMix];
        delete KittyEventOriginal[uMix];
        delete KittyEventRnd[uMix];
        delete KittyEventDimiPhi[uMix];
    }
    delete [] KittyParticleMultRnd;
    delete [] KittyEvent;
    delete [] KittyEventOriginal;
    delete [] KittyEventRnd;
    delete [] KittyEventDimiPhi;
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        delete KittyEventMultRnd[uRnd];
        delete KittyBufferMultRnd[uRnd];
        delete KittyEventDimiPhi_Mult[uRnd];
        delete KittyBufferDimiPhi_Mult[uRnd];
    }
    delete [] KittyBufferOriginal;
    delete [] KittyEventMultRnd;
    delete [] KittyBufferMultRnd;
    delete [] KittyEventDimiPhi_Mult;
    delete [] KittyBufferDimiPhi_Mult;
    delete [] DimiPhi_Diff;

    delete KittyBuffer;
    delete KittyBufferOriginal;
    delete KittyBufferRnd;
    delete KittyBufferDimiPhi;
    delete fOut;
}









void dEta_dPhi_Ck_QS(const TString& TranModDescr, const TString& DataSetDescr, const bool& ALICE_acc){
    const double QS_r = 0.3;
    const double QS_s = 0.4;

    const double kMin = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMax = DataSetDescr=="pp"?6:DataSetDescr=="pLambda"?6:DataSetDescr=="pXim"?6:6;
    const double kMinZoom = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMaxZoom = DataSetDescr=="pp"?2:DataSetDescr=="pLambda"?2:DataSetDescr=="pXim"?2:2;

    const double phiMin = -TMath::Pi()*2;
    const double phiMax = TMath::Pi()*2;

    const double etaMin = -1.5;
    const double etaMax = 1.5;

    const TString OutputFolder = TString::Format("%s/MixedEvents/dEta_dPhi_Ck_QS/",GetFemtoOutputFolder());
    const unsigned NumMomBins = DataSetDescr=="pp"?300:DataSetDescr=="pLambda"?240:DataSetDescr=="pXim"?120:120;
    const unsigned NumPhiBins = NumMomBins;
    const unsigned NumEtaBins = NumMomBins;
    const TString OutFileBaseName = OutputFolder+TranModDescr+"_"+DataSetDescr;
    const TString InputFileName = DataSetDescr=="pp"?TransportFile_pp_Alice:DataSetDescr=="pLambda"?TransportFile_pL_Alice:
        DataSetDescr=="pXim"?TString("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_100KiLrz_pXim.f19"):"";
    const unsigned NumBlankHeaderLines=3;
    const unsigned MaxPairsToRead = 4e6;
    const unsigned HighMultLimit = 128;
    int pdgID[2] = {2212, DataSetDescr=="pp"?2212:DataSetDescr=="pLambda"?3122:DataSetDescr=="pXim"?3312:0};

    const unsigned NumNorm = 1;
    const double kNormMin[NumNorm] = {kMin};
    const double kNormMax[NumNorm] = {kMax};

    TFile* fOut = new TFile(OutFileBaseName+".root","recreate");
/// THE QsRs does not properly work
//to make it work you need to get SE + QS events, than rotate the particles within this particular event
    //hNk in same event
    TH1F* hNkSe = new TH1F("hNkSe", "hNkSe", NumMomBins, kMin, kMax);
    //hNk after including QS
    TH1F* hNkQs = new TH1F("hNkQs", "hNkQs", NumMomBins, kMin, kMax);
    //rotated sample
    TH1F* hNkRs = new TH1F("hNkRs", "hNkRs", NumMomBins, kMin, kMax);
    //rotated sample, after QS
    TH1F* hNkQsRs = new TH1F("hNkQsRs", "hNkQsRs", NumMomBins, kMin, kMax);
    //hNkQs/hNkSe
    TH1F* hCkQs = new TH1F("hCkQs", "hCkQs", NumMomBins, kMin, kMax);
    TH1F* hCkRs = new TH1F("hCkRs", "hCkRs", NumMomBins, kMin, kMax);
    TH1F* hCkQsRs = new TH1F("hCkQsRs", "hCkQsRs", NumMomBins, kMin, kMax);

    //hNk in same event
    TH1F* hNphiSe = new TH1F("hNphiSe", "hNphiSe", NumPhiBins, phiMin, phiMax);
    TH1F* hNphiMe = new TH1F("hNphiMe", "hNphiMe", NumPhiBins, phiMin, phiMax);
    //hNk after including Qs
    TH1F* hNphiQs = new TH1F("hNphiQs", "hNphiQs", NumPhiBins, phiMin, phiMax);
    TH1F* hNphiRs = new TH1F("hNphiRs", "hNphiRs", NumPhiBins, phiMin, phiMax);
    TH1F* hNphiQsRs = new TH1F("hNphiQsRs", "hNphiQsRs", NumPhiBins, phiMin, phiMax);
    //hNkSe/hNkMe
    TH1F* hCphiMe = new TH1F("hCphiMe", "hCphiMe", NumPhiBins, phiMin, phiMax);
    TH1F* hCphiQs = new TH1F("hCphiQs", "hCphiQs", NumPhiBins, phiMin, phiMax);
    TH1F* hCphiRs = new TH1F("hCphiRs", "hCphiRs", NumPhiBins, phiMin, phiMax);
    TH1F* hCphiQsRs = new TH1F("hCphiQsRs", "hCphiQsRs", NumPhiBins, phiMin, phiMax);

    //hNk in same event
    TH1F* hNetaSe = new TH1F("hNetaSe", "hNetaSe", NumEtaBins, etaMin, etaMax);
    //hNk after including Qs
    TH1F* hNetaQs = new TH1F("hNetaQs", "hNetaQs", NumEtaBins, etaMin, etaMax);
    TH1F* hNetaRs = new TH1F("hNetaRs", "hNetaRs", NumPhiBins, phiMin, phiMax);
    TH1F* hNetaQsRs = new TH1F("hNetaQsRs", "hNetaQsRs", NumPhiBins, phiMin, phiMax);
    //hNkSe/hNkMe
    TH1F* hCetaQs = new TH1F("hCetaQs", "hCetaQs", NumEtaBins, etaMin, etaMax);
    TH1F* hCetaRs = new TH1F("hCetaRs", "hCetaRs", NumPhiBins, phiMin, phiMax);
    TH1F* hCetaQsRs = new TH1F("hCetaQsRs", "hCetaQsRs", NumPhiBins, phiMin, phiMax);

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    unsigned NumTotalPairs=0;
    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;
    unsigned uBuffer=0;
    unsigned NumSePairs;
    unsigned NumSePairsRs;

    CatsParticle KittyParticle;
    CatsParticle KittyParticleRs;
    CatsEvent* KittyEvent = new CatsEvent(pdgID[0],pdgID[1]);
    CatsEvent* KittyEventLab = new CatsEvent(pdgID[0],pdgID[1]);
    CatsEvent* KittyEventRs = new CatsEvent(pdgID[0],pdgID[1]);
    CatsEvent* KittyEventRsLab = new CatsEvent(pdgID[0],pdgID[1]);
    CatsDataBuffer* KittyBuffer = new CatsDataBuffer(1,pdgID[0],pdgID[1]);
    CatsDataBuffer* KittyBufferLab = new CatsDataBuffer(1,pdgID[0],pdgID[1]);
    CatsDataBuffer* KittyBufferRs = new CatsDataBuffer(1,pdgID[0],pdgID[1]);
    CatsDataBuffer* KittyBufferRsLab = new CatsDataBuffer(1,pdgID[0],pdgID[1]);

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;
    TRandom3 rangen(11);
    while(!feof(InFile)){
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
        TotNumEvents++;
        if(NumPartInEvent>int(HighMultLimit)) RejectedHighMultEvents++;
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            KittyParticleRs = KittyParticle;
            if(NumTotalPairs>=MaxPairsToRead) continue;
            if(NumPartInEvent>int(HighMultLimit)) continue;

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            //ALICE ACCEPTANCE
            if( (ALICE_acc&&KittyParticle.GetP()>0.4&&fabs(KittyParticle.GetPseudoRap())<1.3) || !ALICE_acc ){
                KittyEvent->AddParticle(KittyParticle);
                KittyEventLab->AddParticle(KittyParticle);
            }
            //for(int iDepth=0; iDepth<8; iDepth++){
            KittyParticleRs.RotateMomPhi(rangen.Uniform(2.*TMath::Pi()));
            if( (ALICE_acc&&KittyParticleRs.GetP()>0.4&&fabs(KittyParticleRs.GetPseudoRap())<1.3) || !ALICE_acc ){
                KittyEventRs->AddParticle(KittyParticleRs);
                KittyEventRsLab->AddParticle(KittyParticleRs);
            }
            //}

        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)

        KittyEvent->ComputeParticlePairs();
        KittyEventLab->ComputeParticlePairs(false,false);
        KittyBuffer->SetEvent(0, *KittyEvent);
        KittyBufferLab->SetEvent(0, *KittyEventLab);

        KittyEventRs->ComputeParticlePairs();
        KittyEventRsLab->ComputeParticlePairs(false,false);
        KittyBufferRs->SetEvent(0, *KittyEventRs);
        KittyBufferRsLab->SetEvent(0, *KittyEventRsLab);

        //if the buffer is full -> empty it!
        //note that if it happens that we leave the while loop before emptying the buffer,
        //uBuffer will be != than zero! use this condition to empty the buffer when exiting the loop!

        KittyBuffer->GoBabyGo();
        KittyBufferLab->GoBabyGo();
        KittyBufferRs->GoBabyGo();
        KittyBufferRsLab->GoBabyGo();
        NumSePairs = KittyBuffer->GetNumPairsSameEvent();
        NumSePairsRs = KittyBufferRs->GetNumPairsSameEvent();
        NumTotalPairs += KittyBuffer->GetNumPairs();
        for(unsigned uSe=0; uSe<NumSePairs; uSe++){
            int PART1 = rangen.Integer(2);
            int PART2 = !PART1;
            hNkSe->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
            hCkRs->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
            hNphiSe->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
            hCphiMe->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
            hCphiRs->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
            hNetaSe->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPseudoRap()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPseudoRap());
            hCetaRs->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPseudoRap()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPseudoRap());
            double Pair_k = KittyBuffer->GetSePair(uSe)->GetP()*0.5;
            double ProbToReject = QS_s*exp(-Pair_k*Pair_k*QS_r*QS_r/0.197327/0.197327);
            if(rangen.Uniform(1)>ProbToReject){
                hNkQs->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkQs->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hNphiQs->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
                hCphiQs->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
                hNetaQs->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPseudoRap()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPseudoRap());
                hCetaQs->Fill(KittyBufferLab->GetSePair(uSe)->GetParticle(PART1).GetPseudoRap()-KittyBufferLab->GetSePair(uSe)->GetParticle(PART2).GetPseudoRap());
            }
        }
        for(unsigned uSe=0; uSe<NumSePairsRs; uSe++){
            int PART1 = rangen.Integer(2);
            int PART2 = !PART1;
            hNkRs->Fill(KittyBufferRs->GetSePair(uSe)->GetP()*0.5);
            hNphiRs->Fill(KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
            hNetaRs->Fill(KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART1).GetPseudoRap()-KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART2).GetPseudoRap());
            double Pair_k = KittyBufferRs->GetSePair(uSe)->GetP()*0.5;
            double ProbToReject = QS_s*exp(-Pair_k*Pair_k*QS_r*QS_r/0.197327/0.197327);
            if(rangen.Uniform(1)>ProbToReject){
                hNkQsRs->Fill(KittyBufferRs->GetSePair(uSe)->GetP()*0.5);
                hCkQsRs->Fill(KittyBufferRs->GetSePair(uSe)->GetP()*0.5);
                hNphiQsRs->Fill(KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
                //hCphiQsRs->Fill(KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART1).GetPphi()-KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART2).GetPphi());
                hNetaQsRs->Fill(KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART1).GetPseudoRap()-KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART2).GetPseudoRap());
                hCetaQsRs->Fill(KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART1).GetPseudoRap()-KittyBufferRsLab->GetSePair(uSe)->GetParticle(PART2).GetPseudoRap());
            }
        }
        for(unsigned uMe=0; uMe<KittyBuffer->GetNumPairsMixedEvent(); uMe++){
            int PART1 = rangen.Integer(2);
            int PART2 = !PART1;
            hNphiMe->Fill(KittyBufferLab->GetMePair(uMe)->GetParticle(PART1).GetPphi()-KittyBufferLab->GetMePair(uMe)->GetParticle(PART2).GetPphi());
        }
        KittyEvent->Reset();
        KittyEventLab->Reset();
        KittyEventRs->Reset();
        KittyEventRsLab->Reset();

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }

    }//while(!feof(InFile))

    if(ProgressBar){
        printf("\r\033[K");
    }


    //for(int iBin=1; iBin<=hCkQs->GetNbinsX(); iBin++){hCkQs->SetBinError(iBin,hNkSe->GetBinError(iBin));}
    //for(int iBin=1; iBin<=hCphiQs->GetNbinsX(); iBin++){hCphiQs->SetBinError(iBin,hNphiSe->GetBinError(iBin));}
    //for(int iBin=1; iBin<=hCetaQs->GetNbinsX(); iBin++){hCetaQs->SetBinError(iBin,hNetaSe->GetBinError(iBin));}

    hNkSe->Sumw2();
    //for(int iBin=1; iBin<=hNkSe->GetNbinsX(); iBin++){hNkSe->SetBinError(iBin,0);}
    hNkQs->Sumw2();
    hNkRs->Sumw2();
    hNkQsRs->Sumw2();
    hCkQs->Sumw2();
    hCkRs->Sumw2();
    hCkQsRs->Sumw2();
    hNphiSe->Sumw2();
    hNphiMe->Sumw2();
    //for(int iBin=1; iBin<=hNphiSe->GetNbinsX(); iBin++){hNphiSe->SetBinError(iBin,0);}
    hNphiQs->Sumw2();
    hNphiRs->Sumw2();
    hNphiQsRs->Sumw2();
    hCphiMe->Sumw2();
    hCphiQs->Sumw2();
    hCphiRs->Sumw2();
    //hCphiQsRs->Sumw2();
    hNetaSe->Sumw2();
    //for(int iBin=1; iBin<=hNetaSe->GetNbinsX(); iBin++){hNetaSe->SetBinError(iBin,0);}
    hNetaQs->Sumw2();
    hNetaRs->Sumw2();
    hNetaQsRs->Sumw2();
    hCetaQs->Sumw2();
    hCetaRs->Sumw2();
    hCetaQsRs->Sumw2();

    hNkSe->Scale(1./hNkSe->Integral(),"width");
    hNkQs->Scale(1./hNkQs->Integral(),"width");
    hNkRs->Scale(1./hNkRs->Integral(),"width");
    hNkQsRs->Scale(1./hNkQsRs->Integral(),"width");
    hCkQs->Scale(1./hCkQs->Integral(),"width");
    hCkRs->Scale(1./hCkRs->Integral(),"width");
    hCkQsRs->Scale(1./hCkQsRs->Integral(),"width");
    hNphiSe->Scale(1./hNphiSe->Integral(),"width");
    hNphiMe->Scale(1./hNphiMe->Integral(),"width");
    hNphiQs->Scale(1./hNphiQs->Integral(),"width");
    hNphiRs->Scale(1./hNphiRs->Integral(),"width");
    hNphiQsRs->Scale(1./hNphiQsRs->Integral(),"width");
    hCphiMe->Scale(1./hCphiMe->Integral(),"width");
    hCphiQs->Scale(1./hCphiQs->Integral(),"width");
    hCphiRs->Scale(1./hCphiRs->Integral(),"width");
    //hCphiQsRs->Scale(1./hCphiQsRs->Integral(),"width");
    hNetaSe->Scale(1./hNetaSe->Integral(),"width");
    hNetaQs->Scale(1./hNetaQs->Integral(),"width");
    hNetaRs->Scale(1./hNetaRs->Integral(),"width");
    hNetaQsRs->Scale(1./hNetaQsRs->Integral(),"width");
    hCetaQs->Scale(1./hCetaQs->Integral(),"width");
    hCetaRs->Scale(1./hCetaRs->Integral(),"width");
    hCetaQsRs->Scale(1./hCetaQsRs->Integral(),"width");

    printf("hNkSe->Integral = %f\n",hNkSe->Integral());
    printf("hNkQs->Integral = %f\n",hNkQs->Integral());
    printf("hCkQs->Integral = %f\n",hCkQs->Integral());

    hCkQs->Divide(hNkSe);
    hCphiMe->Divide(hNphiMe);
    hCphiQs->Divide(hNphiSe);
    hCetaQs->Divide(hNetaSe);

    hCkRs->Divide(hNkRs);
    hCphiRs->Divide(hNphiRs);
    hCetaRs->Divide(hNetaRs);

    hCkQsRs->Divide(hNkRs);
    //hCphiQsRs->Divide(hNphiRs);
    for(int iBin=1; iBin<=hCphiQsRs->GetNbinsX(); iBin++){
        hCphiQsRs->SetBinContent(iBin,hCphiQs->GetBinContent(iBin)*hCphiRs->GetBinContent(iBin));
        hCphiQsRs->SetBinError(iBin,hCphiQs->GetBinError(iBin));
    }

    hCetaQsRs->Divide(hNetaRs);



    TF1* UnitLine = new TF1("UnitLine", "1", kMin, kMax);
    UnitLine->SetLineColor(kGray);
    UnitLine->SetLineWidth(5);

    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumMomBins, kMin, kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (GeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.15);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    //hAxis->GetYaxis()->SetTitle("(N per bin)x10^{-3}");
    hAxis->GetYaxis()->SetTitle("Arbitrary units");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.45);
    hAxis->GetYaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLimits(0, hNkSe->GetBinContent(hNkSe->GetMaximumBin())*1.4);
    hAxis->GetYaxis()->SetRangeUser(0, hNkSe->GetBinContent(hNkSe->GetMaximumBin())*1.4);

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", NumMomBins, kMin, kMax);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio->GetXaxis()->SetTitle("k (GeV)");
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.15);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio->GetYaxis()->SetTitle("N_{SE}/N_{ME}");
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(0.45/1.0);
    hAxisRatio->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
    hAxisRatio->GetYaxis()->SetLimits(0.0, 2.0);
    hAxisRatio->GetYaxis()->SetNdivisions(504);

    TPaveText* InfoText1 = new TPaveText(0.1,0.8,0.7,0.94, "blNDC");//lbrt
    InfoText1->SetName("InfoText1");
    InfoText1->SetBorderSize(1);
    InfoText1->SetTextSize(0.065);
    InfoText1->SetFillColor(kWhite);
    InfoText1->SetTextFont(42);
    InfoText1->AddText("N.B. no multiplicity binning so far!");
    //InfoText1->AddText(TString::Format("Norm. region: [%.2f, %.2f] GeV", kNormMin, kNormMax));



    for(unsigned uNorm=0; uNorm<NumNorm; uNorm++){

        TF1* NormLine = new TF1(TString::Format("NormLine%u",uNorm), "[0]", kNormMin[uNorm], kNormMax[uNorm]);
        NormLine->SetLineColor(kGray);
        NormLine->SetLineWidth(5);
        NormLine->SetParameter(0,1);
        //0 = SE/Qs
        const unsigned NumTypes=1;

        TH1F** HistCk = new TH1F* [NumTypes];
        TH1F** HistSe = new TH1F* [NumTypes];
        TH1F** HistMe = new TH1F* [NumTypes];
        HistCk[0] = hCkQs;
        HistSe[0] = hNkQs;
        HistMe[0] = hNkSe;

        for(unsigned uType=0; uType<NumTypes; uType++){

            //HistCk[uType]->Fit(NormLine, "S, N, R, M");
            //HistSe[uType]->Scale(1./NormLine->GetParameter(0));
            //HistCk[uType]->Scale(1./NormLine->GetParameter(0));

            HistSe[uType]->SetLineWidth(6);
            HistSe[uType]->SetLineColor(kBlue+1);

            HistMe[uType]->SetLineWidth(6);
            HistMe[uType]->SetLineColor(kRed+1);

            HistCk[uType]->SetLineWidth(6);
            HistCk[uType]->SetLineColor(kBlack);

            TH1F* hShade = new TH1F(TString::Format("hShade%u_%u",uNorm,uType), TString::Format("hShade%u_%u",uNorm,uType), 1, kNormMin[uNorm], kNormMax[uNorm]);
            double hShadeMaxVal = HistSe[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.4;
            if(HistCk[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.4>hShadeMaxVal) hShadeMaxVal = HistCk[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.4;
            hShade->SetBinContent(1, hShadeMaxVal);
            hShade->SetFillStyle(3344);//4050
            hShade->SetFillColor(46);
            hShade->SetLineColor(46);
            hShade->SetLineWidth(0);

            TLegend* myLegend = new TLegend(0.7,0.5,0.96,0.94);//lbrt
            myLegend->SetName(TString::Format("myLegend%u_%u",uNorm,uType));
            myLegend->SetTextSize(0.07);
            if(uType==0){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS same event (Qs)");
            }

            myLegend->AddEntry(HistCk[uType], "Ratio");
            myLegend->AddEntry(hShade, "Norm. region");

            DLM_SubPads DrawBoard(1920,1080);
            DrawBoard.AddSubPadTL(0,1,0,0.5);//lrtb
            DrawBoard.AddSubPadTL(0,1,0.5,1);//lrtb

            DrawBoard.SetMargin(0, 0.1, 0.04, 0.0, 0.05);//lrbt
            DrawBoard.SetMargin(1, 0.1, 0.04, 0.15, 0);//lrbt

        ////////////
            DrawBoard.cd(0);

            DrawBoard.SetLabelSize(0, hAxis->GetXaxis(), 18);
            DrawBoard.SetLabelSize(0, hAxis->GetYaxis(), 18);
            DrawBoard.SetTitleSize(0, hAxis->GetXaxis(), 18);
            DrawBoard.SetTitleSize(0, hAxis->GetYaxis(), 18);

            hAxis->GetYaxis()->SetRangeUser(0, HistSe[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.4);
            hAxis->GetXaxis()->SetRangeUser(kMin,kMax);
            hAxisRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
            hAxisRatio->GetXaxis()->SetRangeUser(kMin,kMax);

            hAxis->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            HistSe[uType]->Draw("same");
            HistMe[uType]->Draw("same");
            myLegend->Draw("same");
            InfoText1->Draw("same");

            DrawBoard.cd(1);

            DrawBoard.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoard.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
            DrawBoard.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoard.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

            hAxisRatio->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            UnitLine->Draw("same");
            HistCk[uType]->Draw("same");

            Width_t width_hShade=hShade->GetLineWidth();
            Width_t width_hNkSe=HistSe[uType]->GetLineWidth();
            Width_t width_hNkMe=HistMe[uType]->GetLineWidth();
            Width_t width_UnitLine=UnitLine->GetLineWidth();
            Width_t width_hCk=HistCk[uType]->GetLineWidth();

            DrawBoard.GetCanvas()->SaveAs(TString::Format("%s_DB_N%.0f-%.0f_%u.png",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade/2.5);
            HistSe[uType]->SetLineWidth(width_hNkSe/2.5);
            HistMe[uType]->SetLineWidth(width_hNkMe/2.5);
            UnitLine->SetLineWidth(width_UnitLine/2.5);
            HistCk[uType]->SetLineWidth(width_hCk/2.5);
            DrawBoard.GetCanvas()->SaveAs(TString::Format("%s_DB_N%.0f-%.0f_%u.pdf",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade);
            HistSe[uType]->SetLineWidth(width_hNkSe);
            HistMe[uType]->SetLineWidth(width_hNkMe);
            UnitLine->SetLineWidth(width_UnitLine);
            HistCk[uType]->SetLineWidth(width_hCk);

        /////////////////////////////////////////////////////////////////////
            DLM_SubPads DrawBoardZoom(1920,1080);
            DrawBoardZoom.AddSubPadTL(0,1,0,0.5);//lrtb
            DrawBoardZoom.AddSubPadTL(0,1,0.5,1);//lrtb

            DrawBoardZoom.SetMargin(0, 0.1, 0.04, 0.0, 0.05);//lrbt
            DrawBoardZoom.SetMargin(1, 0.1, 0.04, 0.15, 0);//lrbt

        ////////////
            DrawBoardZoom.cd(0);

            DrawBoardZoom.SetLabelSize(0, hAxis->GetXaxis(), 18);
            DrawBoardZoom.SetLabelSize(0, hAxis->GetYaxis(), 18);
            DrawBoardZoom.SetTitleSize(0, hAxis->GetXaxis(), 18);
            DrawBoardZoom.SetTitleSize(0, hAxis->GetYaxis(), 18);

            unsigned BinLow = HistSe[uType]->FindBin(kMinZoom);
            unsigned BinUp = HistSe[uType]->FindBin(kMaxZoom);
            double MaxYval=0;
            for(unsigned uBin=BinLow; uBin<=BinUp; uBin++){
                if(MaxYval<HistSe[uType]->GetBinContent(uBin)) MaxYval=HistSe[uType]->GetBinContent(uBin);
            }
            hAxis->GetYaxis()->SetRangeUser(0, MaxYval*1.4);
            hAxis->GetYaxis()->SetLimits(0, MaxYval*1.4);
            hAxis->GetXaxis()->SetRangeUser(kMinZoom, kMaxZoom);
            //hAxis->GetXaxis()->SetLimits(kMinZoom, kMaxZoom);

            hAxisRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
            hAxisRatio->GetYaxis()->SetLimits(0.0, 2.0);
            hAxisRatio->GetXaxis()->SetRangeUser(kMinZoom, kMaxZoom);
            //hAxisRatio->GetXaxis()->SetLimits(kMinZoom, kMaxZoom);

            hAxis->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            HistSe[uType]->Draw("same");
            HistMe[uType]->Draw("same");
            myLegend->Draw("same");
            InfoText1->Draw("same");

            DrawBoardZoom.cd(1);

            DrawBoardZoom.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoardZoom.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
            DrawBoardZoom.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoardZoom.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

            hAxisRatio->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            UnitLine->Draw("same");
            HistCk[uType]->Draw("same");

            DrawBoardZoom.GetCanvas()->SaveAs(TString::Format("%s_DBZoom_N%.0f-%.0f_%u.png",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade/2.5);
            HistSe[uType]->SetLineWidth(width_hNkSe/2.5);
            HistMe[uType]->SetLineWidth(width_hNkMe/2.5);
            UnitLine->SetLineWidth(width_UnitLine/2.5);
            HistCk[uType]->SetLineWidth(width_hCk/2.5);
            DrawBoardZoom.GetCanvas()->SaveAs(TString::Format("%s_DBZoom_N%.0f-%.0f_%u.pdf",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade);
            HistSe[uType]->SetLineWidth(width_hNkSe);
            HistMe[uType]->SetLineWidth(width_hNkMe);
            UnitLine->SetLineWidth(width_UnitLine);
            HistCk[uType]->SetLineWidth(width_hCk);

    printf("HistSe[uType]->Integral = %f\n",HistSe[uType]->Integral());
    printf("HistMe[uType]->Integral = %f\n",HistMe[uType]->Integral());
    printf("HistCk[uType]->Integral = %f\n",HistCk[uType]->Integral());

            //HistSe[uType]->Write();
            //HistMe[uType]->Write();
            //HistCk[uType]->Write();

            delete hShade;
            delete myLegend;
        }

        delete [] HistCk;
        delete [] HistSe;
        delete [] HistMe;

        delete NormLine;

    }

    hNkSe->Write();
    hNkQs->Write();
    hNkRs->Write();
    hNkQsRs->Write();
    hCkQs->Write();
    hCkRs->Write();
    hCkQsRs->Write();

    hNphiSe->Write();
    hNphiMe->Write();
    hNphiQs->Write();
    hNphiRs->Write();
    hNphiQsRs->Write();
    hCphiMe->Write();
    hCphiQs->Write();
    hCphiRs->Write();
    hCphiQsRs->Write();

    hNetaSe->Write();
    hNetaQs->Write();
    hNetaRs->Write();
    hNetaQsRs->Write();
    hCetaQs->Write();
    hCetaRs->Write();
    hCetaQsRs->Write();

    delete hNkSe;
    delete hNkQs;
    delete hNkRs;
    delete hNkQsRs;
    delete hCkQs;
    delete hCkRs;
    delete hCkQsRs;

    delete hNphiSe;
    delete hNphiMe;
    delete hNphiQs;
    delete hNphiRs;
    delete hNphiQsRs;
    delete hCphiMe;
    delete hCphiQs;
    delete hCphiRs;
    delete hCphiQsRs;

    delete hNetaSe;
    delete hNetaQs;
    delete hNetaRs;
    delete hNetaQsRs;
    delete hCetaQs;
    delete hCetaRs;
    delete hCetaQsRs;

    delete UnitLine;
    delete hAxis;
    delete hAxisRatio;
    delete InfoText1;
    //delete OutRoot;
    delete [] cdummy;
    delete KittyEvent;
    delete KittyEventLab;
    delete KittyBuffer;
    delete KittyBufferLab;
    delete fOut;
}






/*

void DifferentTechniquesVer1(const TString& DataSetDescr){
    const double kMin = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMax = DataSetDescr=="pp"?3:DataSetDescr=="pLambda"?3:DataSetDescr=="pXim"?3:3;
    const double kMinZoom = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMaxZoom = DataSetDescr=="pp"?0.6:DataSetDescr=="pLambda"?0.6:DataSetDescr=="pXim"?0.6:0.6;
    //const double kNormMin = 0.400;
    //const double kNormMax = 0.600;
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/DifferentTechniquesVer1/";
    const unsigned NumMomBins = DataSetDescr=="pp"?150:DataSetDescr=="pLambda"?75:DataSetDescr=="pXim"?75:75;
    const TString OutFileBaseName = OutputFolder+TranModDescr+"_"+DataSetDescr;
    const TString InputFileName = DataSetDescr=="pp"?TransportFile_pp_Alice:DataSetDescr=="pLambda"?TransportFile_pL_Alice:
        DataSetDescr=="pXim"?TString("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_100KiLrz_pXim.f19"):"";
    const unsigned NumBlankHeaderLines=3;
    const unsigned MaxPairsToRead = 4e7;
    const unsigned HighMultLimit = 128;
    const unsigned MixingDepth = 8;
    const unsigned RandomDepth = 16;
    //0 = ME, 1 = DimiPhi, 2 = RandomPhi, 3 = Stavinski
    const unsigned NumMethods = 3;
    TString MethodDescription[NumMethods] = {"ME,DimiPhi,RandPhi"};
    const double DimiSmearPhi = (10.)*DegToRad;
    int pdgID[2] = {2212, DataSetDescr=="pp"?2212:DataSetDescr=="pLambda"?3122:DataSetDescr=="pXim"?3312:0};

    //hNk in same event
    TH1F* hNkSe = new TH1F("hNkSe", "hNkSe", NumMomBins, kMin, kMax);
    //hNk in the reference of each model
    TH1F** hNkRe = new TH1F* [NumMethods];
    TH1F** hCk = new TH1F* [NumMethods];
    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        hNkRe[uMeth] = new TH1F("hNkRe_"+MethodDescription[uMeth], "hNkRe_"+MethodDescription[uMeth], NumMomBins, kMin, kMax);
        hCk[uMeth] = new TH1F("hCk_"+MethodDescription[uMeth], "hCk_"+MethodDescription[uMeth], NumMomBins, kMin, kMax);
    }

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    unsigned NumTotalPairs=0;
    //unsigned NumTotalPairsRnd=0;
    unsigned NumAllPairsRnd=0;
    unsigned NumAllPairsDimiPhi=0;
    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;
    unsigned uBuffer=0;
    unsigned NumSePairs;
    unsigned NumMePairs;
    unsigned NumSePairsRnd;
    unsigned NumSePairsDimiPhi;
    //unsigned NumMePairsRnd;
    unsigned NumSePairsMultRnd;
    unsigned NumSePairsDimiPhi_Mult;

    CatsParticle KittyParticle;
    CatsParticle KittyParticleDimiPhi;

    CatsParticle* KittyParticleMultRnd = new CatsParticle[RandomDepth];

    CatsParticle KittyParticleStavinsky;
    //rotated in the phi plane (only a bit). The next particle -> rotation fixed such as to conserve E-P.
    CatsParticle KittyParticleDimiPhi;
    //same but repeated multiple times
    CatsParticle* KittyParticleDimiPhi_Mult = new CatsParticle[RandomDepth];

    CatsEvent** KittyEvent = new CatsEvent* [MixingDepth];
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        KittyEvent[uMix] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    CatsEvent** KittyEventRnd = new CatsEvent* [MixingDepth];
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        KittyEventRnd[uMix] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    //this has different structure on purpose
    CatsEvent** KittyEventMultRnd = new CatsEvent* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyEventMultRnd[uRnd] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    CatsEvent** KittyEventDimiPhi = new CatsEvent* [MixingDepth];
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        KittyEventDimiPhi[uMix] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    //this has different structure on purpose
    CatsEvent** KittyEventDimiPhi_Mult = new CatsEvent* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyEventDimiPhi_Mult[uRnd] = new CatsEvent(pdgID[0],pdgID[1]);
    }



    CatsDataBuffer* KittyBuffer = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    CatsDataBuffer* KittyBufferRnd = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    CatsDataBuffer** KittyBufferMultRnd = new CatsDataBuffer* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyBufferMultRnd[uRnd] = new CatsDataBuffer(1,pdgID[0],pdgID[1]);
    }
    CatsDataBuffer* KittyBufferDimiPhi = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    CatsDataBuffer** KittyBufferDimiPhi_Mult = new CatsDataBuffer* [RandomDepth];
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        KittyBufferDimiPhi_Mult[uRnd] = new CatsDataBuffer(1,pdgID[0],pdgID[1]);
    }
    CatsLorentzVector* DimiPhi_Diff = new CatsLorentzVector[RandomDepth];

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;
    TRandom3 rangen(11);
    while(!feof(InFile)){
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
        TotNumEvents++;
        if(NumPartInEvent>int(HighMultLimit)) RejectedHighMultEvents++;

        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            KittyParticleRnd = KittyParticle;
            KittyParticleRnd.RotateMomPhi(rangen.Uniform(2.*TMath::Pi()));

            //KittyParticleDimi1 = KittyParticle
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                KittyParticleMultRnd[uRnd] = KittyParticle;
                KittyParticleMultRnd[uRnd].RotateMomPhi(rangen.Uniform(2.*TMath::Pi()));

                KittyParticleDimiPhi_Mult[uRnd] = KittyParticle;
            }
            if(NumTotalPairs>=MaxPairsToRead) continue;
            if(NumPartInEvent>int(HighMultLimit)) continue;

            //ALICE ACCEPTANCE
            if(KittyParticle.GetP()<0.4) continue;
            if(fabs(KittyParticle.GetPseudoRap())>0.8) continue;//!

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            KittyEvent[uBuffer]->AddParticle(KittyParticle);
            KittyEventRnd[uBuffer]->AddParticle(KittyParticleRnd);
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                KittyEventMultRnd[uRnd]->AddParticle(KittyParticleMultRnd[uRnd]);
            }

            //Dimi's random phi
            //rotation are done pair-wise, such that each next rotation restores E-Mom conservation
            //the last particle is not rotated, as there are no particles left to compensate
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                if(iPart%2==0&&iPart!=NumPartInEvent){
                    KittyParticleDimiPhi_Mult[uRnd].RotateMomPhi(rangen.Uniform(DimiSmearPhi));
                    DimiPhi_Diff[uRnd] = KittyParticleDimiPhi_Mult[uRnd]-KittyParticle;
                }
                else if(iPart%2==1){
                    //diff = new-old
                    //new = diff+old
                    //=> new = old-diff (we want to change with -diff)
                    KittyParticleDimiPhi_Mult[uRnd] = KittyParticle-DimiPhi_Diff[uRnd];
                }
                else{
                    //nothing
                }
            }
            KittyEventDimiPhi[uBuffer]->AddParticle(KittyParticleDimiPhi_Mult[0]);
            for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
                KittyEventDimiPhi_Mult[uRnd]->AddParticle(KittyParticleDimiPhi_Mult[uRnd]);
            }
        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)

        KittyEvent[uBuffer]->ComputeParticlePairs();
        KittyBuffer->SetEvent(uBuffer, *KittyEvent[uBuffer]);

        KittyEventRnd[uBuffer]->ComputeParticlePairs();
        KittyBufferRnd->SetEvent(uBuffer, *KittyEventRnd[uBuffer]);

        KittyEventDimiPhi[uBuffer]->ComputeParticlePairs();
        KittyBufferDimiPhi->SetEvent(uBuffer, *KittyEventDimiPhi[uBuffer]);

        for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
            KittyEventMultRnd[uRnd]->ComputeParticlePairs();
            KittyBufferMultRnd[uRnd]->SetEvent(uBuffer, *KittyEventMultRnd[uBuffer]);

            KittyEventDimiPhi_Mult[uRnd]->ComputeParticlePairs();
            KittyBufferDimiPhi_Mult[uRnd]->SetEvent(uBuffer, *KittyEventDimiPhi_Mult[uBuffer]);
        }

        uBuffer++;

        //if the buffer is full -> empty it!
        //note that if it happens that we leave the while loop before emptying the buffer,
        //uBuffer will be != than zero! use this condition to empty the buffer when exiting the loop!
        if(uBuffer==MixingDepth){
            KittyBuffer->GoBabyGo();
            NumSePairs = KittyBuffer->GetNumPairsSameEvent();
            NumMePairs = KittyBuffer->GetNumPairsMixedEvent();
            NumTotalPairs += KittyBuffer->GetNumPairs();
            for(unsigned uSe=0; uSe<NumSePairs; uSe++){
                hNkSe->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCk->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkRnd->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkMeRnd->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkDimiPhi->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
                hCkDimiPhi_Mult->Fill(KittyBuffer->GetSePair(uSe)->GetP()*0.5);
            }
            for(unsigned uMe=0; uMe<NumMePairs; uMe++){
                hNkMe->Fill(KittyBuffer->GetMePair(uMe)->GetP()*0.5);
                hCkMeMeRnd->Fill(KittyBuffer->GetMePair(uMe)->GetP()*0.5);
            }

            KittyBufferRnd->GoBabyGo();
            NumSePairsRnd = KittyBufferRnd->GetNumPairsSameEvent();
            //NumMePairsRnd = KittyBufferRnd->GetNumPairsMixedEvent();
            NumAllPairsRnd = KittyBufferRnd->GetNumPairs();
            //double AvgNumPairsRnd = KittyBufferRnd->GetAvgNumPairs();
            for(unsigned uSe=0; uSe<NumSePairsRnd; uSe++){
                hNkSeRnd->Fill(KittyBufferRnd->GetSePair(uSe)->GetP()*0.5);
                hCkSeRndMe->Fill(KittyBufferRnd->GetSePair(uSe)->GetP()*0.5);
                hCkSeRndMeRnd->Fill(KittyBufferRnd->GetSePair(uSe)->GetP()*0.5);
            }
            //note that for the Rnd ME we can take ALL pairs, including from SE, since they should all
            //by construction carry no correlations!
            for(unsigned uMe=0; uMe<NumAllPairsRnd; uMe++){
                hNkMeRnd->Fill(KittyBufferRnd->GetPair(uMe)->GetP()*0.5);
            }

            KittyBufferDimiPhi->GoBabyGo();
            NumSePairsDimiPhi = KittyBufferDimiPhi->GetNumPairsSameEvent();
            //NumMePairsDimiPhi = KittyBufferDimiPhi->GetNumPairsMixedEvent();
//printf("NumSePairsDimiPhi=%u\n",NumSePairsDimiPhi);
            NumAllPairsDimiPhi = KittyBufferDimiPhi->GetNumPairs();
            //double AvgNumPairsDimiPhi = KittyBufferDimiPhi->GetAvgNumPairs();
            for(unsigned uSe=0; uSe<NumSePairsDimiPhi; uSe++){
                hNkSeDimiPhi->Fill(KittyBufferDimiPhi->GetSePair(uSe)->GetP()*0.5);
                //hCkSeDimiPhiMe->Fill(KittyBufferDimiPhi->GetSePair(uSe)->GetP()*0.5);
                //hCkSeDimiPhiMeDimiPhi->Fill(KittyBufferDimiPhi->GetSePair(uSe)->GetP()*0.5);
            }
            //note that for the DimiPhi ME we can take ALL pairs, including from SE, since they should all
            //by construction carry no correlations!
            for(unsigned uMe=0; uMe<NumAllPairsDimiPhi; uMe++){
                hNkMeDimiPhi->Fill(KittyBufferDimiPhi->GetPair(uMe)->GetP()*0.5);
            }

            for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
                KittyEvent[uDepth]->Reset();
                KittyEventRnd[uDepth]->Reset();
                KittyEventDimiPhi[uDepth]->Reset();
            }
            uBuffer=0;
        }

        for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
            KittyBufferMultRnd[uRnd]->GoBabyGo();
            NumSePairsMultRnd = KittyBufferMultRnd[uRnd]->GetNumPairsSameEvent();
            for(unsigned uSe=0; uSe<NumSePairsMultRnd; uSe++){
                hNkMultSeRnd->Fill(KittyBufferMultRnd[uRnd]->GetSePair(uSe)->GetP()*0.5);
            }
            KittyEventMultRnd[uRnd]->Reset();

            KittyBufferDimiPhi_Mult[uRnd]->GoBabyGo();
            NumSePairsDimiPhi_Mult = KittyBufferDimiPhi_Mult[uRnd]->GetNumPairsSameEvent();
//printf("NumSePairsDimiPhi_Mult=%u\n",NumSePairsDimiPhi_Mult);
            for(unsigned uSe=0; uSe<NumSePairsDimiPhi_Mult; uSe++){
                hNkSeDimiPhi_Mult->Fill(KittyBufferDimiPhi_Mult[uRnd]->GetSePair(uSe)->GetP()*0.5);
            }
            KittyEventDimiPhi_Mult[uRnd]->Reset();
        }

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }

    }//while(!feof(InFile))




    if(ProgressBar){
        printf("\r\033[K");
    }

    hNkSe->Sumw2();
    hNkMe->Sumw2();
    hCk->Sumw2();

    hNkSeRnd->Sumw2();
    hNkMeRnd->Sumw2();
    hCkRnd->Sumw2();
    hCkMeRnd->Sumw2();
    hCkMeMeRnd->Sumw2();
    hCkSeRndMe->Sumw2();
    hCkSeRndMeRnd->Sumw2();

    hNkMultSeRnd->Sumw2();

    hNkSeDimiPhi->Sumw2();
    hNkSeDimiPhi_Mult->Sumw2();
    hNkMeDimiPhi->Sumw2();
    hCkDimiPhi->Sumw2();
    hCkDimiPhi_Mult->Sumw2();

    TH1F* h_KinColl_Cons = new TH1F("h_KinColl_Cons", "h_KinColl_Cons", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        h_KinColl_Cons->SetBinContent(uBin+1, hNkMe->GetBinContent(uBin+1));
        h_KinColl_Cons->SetBinError(uBin+1, hNkMe->GetBinError(uBin+1));
    }
    h_KinColl_Cons->Divide(hNkMeRnd);
    h_KinColl_Cons->Multiply(hNkSeRnd);
    TH1F* hCk_KinColl_Cons = new TH1F("hCk_KinColl_Cons", "hCk_KinColl_Cons", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        hCk_KinColl_Cons->SetBinContent(uBin+1, hNkSe->GetBinContent(uBin+1));
        hCk_KinColl_Cons->SetBinError(uBin+1, hNkSe->GetBinError(uBin+1));
    }

    //hNkSe/hNkMultSeRnd
    TH1F* hCkMultRnd = new TH1F("hCkMultRnd", "hCkMultRnd", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        hCkMultRnd->SetBinContent(uBin+1, hNkSe->GetBinContent(uBin+1));
        hCkMultRnd->SetBinError(uBin+1, hNkSe->GetBinError(uBin+1));
    }

    //hNkSeRnd/hNkMultSeRnd
    TH1F* hCkMultRnd_SeRnd = new TH1F("hCkMultRnd_SeRnd", "hCkMultRnd_SeRnd", NumMomBins, kMin, kMax);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        hCkMultRnd_SeRnd->SetBinContent(uBin+1, hNkSeRnd->GetBinContent(uBin+1));
        hCkMultRnd_SeRnd->SetBinError(uBin+1, hNkSeRnd->GetBinError(uBin+1));
    }

printf("hNkSe->Integral()=%e\n",hNkSe->Integral());
printf("hNkSeRnd->Integral()=%e\n",hNkSeRnd->Integral());
printf("hNkMultSeRnd->Integral()=%e\n",hNkMultSeRnd->Integral());
printf("hNkSeDimiPhi->Integral()=%e\n",hNkSeDimiPhi->Integral());
printf("hNkSeDimiPhi_Mult->Integral()=%e\n",hNkSeDimiPhi_Mult->Integral());

    hNkSe->Scale(100./hNkSe->Integral());
    hNkMe->Scale(100./hNkMe->Integral());
    hCk->Scale(100./hCk->Integral());

    hNkSeRnd->Scale(100./hNkSeRnd->Integral());
    hNkMeRnd->Scale(100./hNkMeRnd->Integral());
    hCkRnd->Scale(100./hCkRnd->Integral());
    hCkMeRnd->Scale(100./hCkMeRnd->Integral());
    hCkMeMeRnd->Scale(100./hCkMeMeRnd->Integral());
    hCkSeRndMe->Scale(100./hCkSeRndMe->Integral());
    hCkSeRndMeRnd->Scale(100./hCkSeRndMeRnd->Integral());
    hNkMeDimiPhi->Scale(100./hNkMeDimiPhi->Integral());
    hNkSeDimiPhi->Scale(100./hNkSeDimiPhi->Integral());
    hNkSeDimiPhi_Mult->Scale(100./hNkSeDimiPhi_Mult->Integral());
    hCkDimiPhi->Scale(100./hCkDimiPhi->Integral());
    hCkDimiPhi_Mult->Scale(100./hCkDimiPhi_Mult->Integral());

    h_KinColl_Cons->Scale(100./h_KinColl_Cons->Integral());
    hCk_KinColl_Cons->Scale(100./hCk_KinColl_Cons->Integral());

    hNkMultSeRnd->Scale(100./hNkMultSeRnd->Integral());
    hCkMultRnd->Scale(100./hCkMultRnd->Integral());
    hCkMultRnd_SeRnd->Scale(100./hCkMultRnd_SeRnd->Integral());

    hCk->Divide(hNkMe);
    hCkRnd->Divide(hNkSeRnd);
    hCkMeRnd->Divide(hNkMeRnd);
    hCkMeMeRnd->Divide(hNkMeRnd);
    hCkSeRndMe->Divide(hNkMe);
    hCkSeRndMeRnd->Divide(hNkMeRnd);
    hCkMultRnd->Divide(hNkMultSeRnd);
    hCkMultRnd_SeRnd->Divide(hNkMultSeRnd);

    hCkDimiPhi->Divide(hNkSeDimiPhi);
    hCkDimiPhi_Mult->Divide(hNkSeDimiPhi_Mult);

    hCk_KinColl_Cons->Divide(h_KinColl_Cons);

    //TFile* OutRoot = new TFile(TString::Format("%s.root",OutFileBaseName.Data()),"recreate");
    //OutRoot->cd();
    //hNkSe->Write();
    //hNkMe->Write();
    //hCk->Write();

    TF1* UnitLine = new TF1("UnitLine", "1", kMin, kMax);
    UnitLine->SetLineColor(kGray);
    UnitLine->SetLineWidth(5);

    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumMomBins, kMin, kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (GeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.15);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    //hAxis->GetYaxis()->SetTitle("(N per bin)x10^{-3}");
    hAxis->GetYaxis()->SetTitle("Arbitrary units");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.45);
    hAxis->GetYaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLimits(0, hNkSe->GetBinContent(hNkSe->GetMaximumBin())*1.2);
    hAxis->GetYaxis()->SetRangeUser(0, hNkSe->GetBinContent(hNkSe->GetMaximumBin())*1.2);

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", NumMomBins, kMin, kMax);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio->GetXaxis()->SetTitle("k (GeV)");
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.15);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio->GetYaxis()->SetTitle("N_{SE}/N_{ME}");
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(0.45/1.0);
    hAxisRatio->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    hAxisRatio->GetYaxis()->SetLimits(0.8, 1.2);
    hAxisRatio->GetYaxis()->SetNdivisions(504);

    TPaveText* InfoText1 = new TPaveText(0.1,0.8,0.7,0.94, "blNDC");//lbrt
    InfoText1->SetName("InfoText1");
    InfoText1->SetBorderSize(1);
    InfoText1->SetTextSize(0.065);
    InfoText1->SetFillColor(kWhite);
    InfoText1->SetTextFont(42);
    InfoText1->AddText("N.B. no multiplicity binning so far!");
    //InfoText1->AddText(TString::Format("Norm. region: [%.2f, %.2f] GeV", kNormMin, kNormMax));

    TFile* fOut = new TFile(OutFileBaseName+".root","recreate");

    for(unsigned uNorm=0; uNorm<NumNorm; uNorm++){

        TF1* NormLine = new TF1(TString::Format("NormLine%u",uNorm), "[0]", kNormMin[uNorm], kNormMax[uNorm]);
        NormLine->SetLineColor(kGray);
        NormLine->SetLineWidth(5);
        NormLine->SetParameter(0,1);

        //hAxis->GetYaxis()->SetRangeUser(0.75, 1.3);
        //hAxis->GetYaxis()->SetLimits(0.75, 1.3);


        //0 = SE/ME
        //1 = SE/RandSE
        //2 = SE/RandME
        //3 = ME/RandME
        //4 = RandSE/ME
        //5 = RandSE/RandME
        //6 = SE/SomethingThatShouldHaveOnlyCollisionKinematicsAndConservationLaws
        //7 = SE/MultSeRnd
        //8 = RandSE/MultSeRnd
        //9 = SE/DimiPhi
        //10 = SE/DimiPhi_Mult
        const unsigned NumTypes=11;

        TH1F** HistCk = new TH1F* [NumTypes];
        TH1F** HistSe = new TH1F* [NumTypes];
        TH1F** HistMe = new TH1F* [NumTypes];
        HistCk[0] = hCk;
        HistSe[0] = hNkSe;
        HistMe[0] = hNkMe;
        HistCk[1] = hCkRnd;
        HistSe[1] = hNkSe;
        HistMe[1] = hNkSeRnd;
        HistCk[2] = hCkMeRnd;
        HistSe[2] = hNkSe;
        HistMe[2] = hNkMeRnd;
        HistCk[3] = hCkMeMeRnd;
        HistSe[3] = hNkMe;
        HistMe[3] = hNkMeRnd;
        HistCk[4] = hCkSeRndMe;
        HistSe[4] = hNkSeRnd;
        HistMe[4] = hNkMe;
        HistCk[5] = hCkSeRndMeRnd;
        HistSe[5] = hNkSeRnd;
        HistMe[5] = hNkMeRnd;
        HistCk[6] = hCk_KinColl_Cons;
        HistSe[6] = hNkSe;
        HistMe[6] = h_KinColl_Cons;
        HistCk[7] = hCkMultRnd;
        HistSe[7] = hNkSe;
        HistMe[7] = hNkMultSeRnd;
        HistCk[8] = hCkMultRnd_SeRnd;
        HistSe[8] = hNkSeRnd;
        HistMe[8] = hNkMultSeRnd;
        HistCk[9] = hCkDimiPhi;
        HistSe[9] = hNkSe;
        HistMe[9] = hNkSeDimiPhi;
        HistCk[10] = hCkDimiPhi_Mult;
        HistSe[10] = hNkSe;
        HistMe[10] = hNkSeDimiPhi_Mult;

        for(unsigned uType=0; uType<NumTypes; uType++){

            HistCk[uType]->Fit(NormLine, "S, N, R, M");
            HistSe[uType]->Scale(1./NormLine->GetParameter(0));
            HistCk[uType]->Scale(1./NormLine->GetParameter(0));

            HistSe[uType]->SetLineWidth(6);
            HistSe[uType]->SetLineColor(kBlue+1);

            HistMe[uType]->SetLineWidth(6);
            HistMe[uType]->SetLineColor(kRed+1);

            HistCk[uType]->SetLineWidth(6);
            HistCk[uType]->SetLineColor(kBlack);

            TH1F* hShade = new TH1F(TString::Format("hShade%u_%u",uNorm,uType), TString::Format("hShade%u_%u",uNorm,uType), 1, kNormMin[uNorm], kNormMax[uNorm]);
            double hShadeMaxVal = HistSe[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2;
            if(HistCk[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2>hShadeMaxVal) hShadeMaxVal = HistCk[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2;
            hShade->SetBinContent(1, hShadeMaxVal);
            hShade->SetFillStyle(3344);//4050
            hShade->SetFillColor(46);
            hShade->SetLineColor(46);
            hShade->SetLineWidth(0);

            TLegend* myLegend = new TLegend(0.7,0.5,0.96,0.94);//lbrt
            myLegend->SetName(TString::Format("myLegend%u_%u",uNorm,uType));
            myLegend->SetTextSize(0.07);
            if(uType==0){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS mixed event");
            }
            else if(uType==1){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized SE");
            }
            else if(uType==2){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized ME");
            }
            else if(uType==3){
                myLegend->AddEntry(HistSe[uType], "EPOS mixed event");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized ME");
            }
            else if(uType==4){
                myLegend->AddEntry(HistSe[uType], "EPOS randomized SE");
                myLegend->AddEntry(HistMe[uType], "EPOS ME");
            }
            else if(uType==5){
                myLegend->AddEntry(HistSe[uType], "EPOS randomized SE");
                myLegend->AddEntry(HistMe[uType], "EPOS randomized ME");
            }
            else if(uType==6){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Dimi's fantasy");
            }
            else if(uType==7){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Andi's fantasy");
            }
            else if(uType==8){
                myLegend->AddEntry(HistSe[uType], "EPOS randomized SE");
                myLegend->AddEntry(HistMe[uType], "Andi's fantasy");
            }
            else if(uType==9){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Dimi random phi");
            }
            else if(uType==10){
                myLegend->AddEntry(HistSe[uType], "EPOS same event");
                myLegend->AddEntry(HistMe[uType], "Dimi random phi (mult)");
            }

            myLegend->AddEntry(HistCk[uType], "Ratio");
            myLegend->AddEntry(hShade, "Norm. region");

            DLM_SubPads DrawBoard(1920,1080);
            DrawBoard.AddSubPadTL(0,1,0,0.5);//lrtb
            DrawBoard.AddSubPadTL(0,1,0.5,1);//lrtb

            DrawBoard.SetMargin(0, 0.1, 0.04, 0.0, 0.05);//lrbt
            DrawBoard.SetMargin(1, 0.1, 0.04, 0.15, 0);//lrbt

        ////////////
            DrawBoard.cd(0);

            DrawBoard.SetLabelSize(0, hAxis->GetXaxis(), 18);
            DrawBoard.SetLabelSize(0, hAxis->GetYaxis(), 18);
            DrawBoard.SetTitleSize(0, hAxis->GetXaxis(), 18);
            DrawBoard.SetTitleSize(0, hAxis->GetYaxis(), 18);

            hAxis->GetYaxis()->SetRangeUser(0, HistSe[uType]->GetBinContent(HistSe[uType]->GetMaximumBin())*1.2);
            hAxis->GetXaxis()->SetRangeUser(kMin,kMax);
            hAxisRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
            hAxisRatio->GetXaxis()->SetRangeUser(kMin,kMax);

            hAxis->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            HistSe[uType]->Draw("same");
            HistMe[uType]->Draw("same");
            myLegend->Draw("same");
            InfoText1->Draw("same");

            DrawBoard.cd(1);

            DrawBoard.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoard.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
            DrawBoard.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoard.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

            hAxisRatio->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            UnitLine->Draw("same");
            HistCk[uType]->Draw("same");

            Width_t width_hShade=hShade->GetLineWidth();
            Width_t width_hNkSe=HistSe[uType]->GetLineWidth();
            Width_t width_hNkMe=HistMe[uType]->GetLineWidth();
            Width_t width_UnitLine=UnitLine->GetLineWidth();
            Width_t width_hCk=HistCk[uType]->GetLineWidth();

            DrawBoard.GetCanvas()->SaveAs(TString::Format("%s_DB_N%.0f-%.0f_%u.png",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade/2.5);
            HistSe[uType]->SetLineWidth(width_hNkSe/2.5);
            HistMe[uType]->SetLineWidth(width_hNkMe/2.5);
            UnitLine->SetLineWidth(width_UnitLine/2.5);
            HistCk[uType]->SetLineWidth(width_hCk/2.5);
            DrawBoard.GetCanvas()->SaveAs(TString::Format("%s_DB_N%.0f-%.0f_%u.pdf",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade);
            HistSe[uType]->SetLineWidth(width_hNkSe);
            HistMe[uType]->SetLineWidth(width_hNkMe);
            UnitLine->SetLineWidth(width_UnitLine);
            HistCk[uType]->SetLineWidth(width_hCk);

        /////////////////////////////////////////////////////////////////////
            DLM_SubPads DrawBoardZoom(1920,1080);
            DrawBoardZoom.AddSubPadTL(0,1,0,0.5);//lrtb
            DrawBoardZoom.AddSubPadTL(0,1,0.5,1);//lrtb

            DrawBoardZoom.SetMargin(0, 0.1, 0.04, 0.0, 0.05);//lrbt
            DrawBoardZoom.SetMargin(1, 0.1, 0.04, 0.15, 0);//lrbt

        ////////////
            DrawBoardZoom.cd(0);

            DrawBoardZoom.SetLabelSize(0, hAxis->GetXaxis(), 18);
            DrawBoardZoom.SetLabelSize(0, hAxis->GetYaxis(), 18);
            DrawBoardZoom.SetTitleSize(0, hAxis->GetXaxis(), 18);
            DrawBoardZoom.SetTitleSize(0, hAxis->GetYaxis(), 18);

            unsigned BinLow = HistSe[uType]->FindBin(kMinZoom);
            unsigned BinUp = HistSe[uType]->FindBin(kMaxZoom);
            double MaxYval=0;
            for(unsigned uBin=BinLow; uBin<=BinUp; uBin++){
                if(MaxYval<HistSe[uType]->GetBinContent(uBin)) MaxYval=HistSe[uType]->GetBinContent(uBin);
            }
            hAxis->GetYaxis()->SetRangeUser(0, MaxYval*1.2);
            hAxis->GetYaxis()->SetLimits(0, MaxYval*1.2);
            hAxis->GetXaxis()->SetRangeUser(kMinZoom, kMaxZoom);
            //hAxis->GetXaxis()->SetLimits(kMinZoom, kMaxZoom);

            hAxisRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
            hAxisRatio->GetYaxis()->SetLimits(0.8, 1.2);
            hAxisRatio->GetXaxis()->SetRangeUser(kMinZoom, kMaxZoom);
            //hAxisRatio->GetXaxis()->SetLimits(kMinZoom, kMaxZoom);

            hAxis->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            HistSe[uType]->Draw("same");
            HistMe[uType]->Draw("same");
            myLegend->Draw("same");
            InfoText1->Draw("same");

            DrawBoardZoom.cd(1);

            DrawBoardZoom.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoardZoom.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
            DrawBoardZoom.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
            DrawBoardZoom.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

            hAxisRatio->Draw("AXIS");
            if(kNormMin[uNorm]!=kMin || kNormMax[uNorm]!=kMax) hShade->Draw("same,HIST");
            UnitLine->Draw("same");
            HistCk[uType]->Draw("same");

            DrawBoardZoom.GetCanvas()->SaveAs(TString::Format("%s_DBZoom_N%.0f-%.0f_%u.png",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade/2.5);
            HistSe[uType]->SetLineWidth(width_hNkSe/2.5);
            HistMe[uType]->SetLineWidth(width_hNkMe/2.5);
            UnitLine->SetLineWidth(width_UnitLine/2.5);
            HistCk[uType]->SetLineWidth(width_hCk/2.5);
            DrawBoardZoom.GetCanvas()->SaveAs(TString::Format("%s_DBZoom_N%.0f-%.0f_%u.pdf",OutFileBaseName.Data(),kNormMin[uNorm]*1000,kNormMax[uNorm]*1000,uType));
            hShade->SetLineWidth(width_hShade);
            HistSe[uType]->SetLineWidth(width_hNkSe);
            HistMe[uType]->SetLineWidth(width_hNkMe);
            UnitLine->SetLineWidth(width_UnitLine);
            HistCk[uType]->SetLineWidth(width_hCk);

            HistSe[uType]->Write();
            HistMe[uType]->Write();
            HistCk[uType]->Write();

            delete hShade;
            delete myLegend;
        }

        delete [] HistCk;
        delete [] HistSe;
        delete [] HistMe;

        delete NormLine;

    }

    delete hNkSe;
    delete hNkSeRnd;
    delete hNkMeRnd;
    delete hNkMe;
    delete hNkMeDimiPhi;
    delete hNkSeDimiPhi;
    delete hNkSeDimiPhi_Mult;
    delete hCk;
    delete hCkRnd;
    delete hCkMeRnd;
    delete hCkMeMeRnd;
    delete hCkSeRndMe;
    delete hCkSeRndMeRnd;
    delete hCkMultRnd;
    delete hCkMultRnd_SeRnd;
    delete hCkDimiPhi;
    delete hCkDimiPhi_Mult;

    delete h_KinColl_Cons;
    delete hCk_KinColl_Cons;

    delete UnitLine;
    delete hAxis;
    delete hAxisRatio;
    delete InfoText1;
    //delete OutRoot;

    delete [] cdummy;
    for(unsigned uMix=0; uMix<MixingDepth; uMix++){
        delete KittyEvent[uMix];
        delete KittyEventRnd[uMix];
        delete KittyEventDimiPhi[uMix];
    }
    delete [] KittyParticleMultRnd;
    delete [] KittyEvent;
    delete [] KittyEventRnd;
    delete [] KittyEventDimiPhi;
    for(unsigned uRnd=0; uRnd<RandomDepth; uRnd++){
        delete KittyEventMultRnd[uRnd];
        delete KittyBufferMultRnd[uRnd];
        delete KittyEventDimiPhi_Mult[uRnd];
        delete KittyBufferDimiPhi_Mult[uRnd];
    }
    delete [] KittyEventMultRnd;
    delete [] KittyBufferMultRnd;
    delete [] KittyEventDimiPhi_Mult;
    delete [] KittyBufferDimiPhi_Mult;
    delete [] DimiPhi_Diff;

    delete KittyBuffer;
    delete KittyBufferRnd;
    delete KittyBufferDimiPhi;
    delete fOut;
}

*/




DLM_CkDecomposition* FIT_CK;
//[0] = Radius
//[1] = Stability
//[2]... the spline pars
double Fit_SE_Splines(double* xVal, double* pars){
    double& Momentum = *xVal;
    FIT_CK->GetCk()->SetSourcePar(0,pars[0]);
    if(FIT_CK->GetCk()->GetNumSourcePar()>1) FIT_CK->GetCk()->SetSourcePar(1,pars[1]);
    double CkVal = FIT_CK->EvalCk(Momentum);
    double BlVal = DLM_FITTER2_FUNCTION_SPLINE3(xVal,&pars[2]);
    double FitVal = CkVal*BlVal;
    return FitVal;
}
//[2] norm
double Fit_Ck_Norm(double* xVal, double* pars){
    double& Momentum = *xVal;
    FIT_CK->GetCk()->SetSourcePar(0,pars[0]);
    if(FIT_CK->GetCk()->GetNumSourcePar()>1) FIT_CK->GetCk()->SetSourcePar(1,pars[1]);
    double CkVal = FIT_CK->EvalCk(Momentum);
    double BlVal = pars[2];
    double FitVal = CkVal*BlVal;
    return FitVal;
}

DLM_Histo<double>* HME_FIT;
double Fit_StretchME(double* xVal, double* pars){
    double EVAL = *xVal/pars[0];
    return HME_FIT->Eval(&EVAL)*pars[0];
}

void CompareSameMixedEventToBoltzmann(){
    //The 1D Boltzmann is actually Gauss with var = s^2 = kT/m
    //=> the distribution of relative momenta between particles following a Boltzmann
    //have the same PDF as S(r), only that we have k as an argument! Isn't that Great!
    //=> I should be able to fit the N(k) with the function we use for the Gaussian source!
    // and if it does not fit, why not use a Levy :D

    const TString FittingMode_pp = "Spline3_5";
    const TString InputFileName = "/home/dmihaylov/CernBox/HM13TeV/AnalysisData/ClosePairRej/SelectedPairs/AnalysisResults.root";
    const TString OutputFolderName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/CompareSameMixedEventToBoltzmann/";;

    const double kMin = 0;
    const double kMax = 540;

    TFile* fInput = new TFile(InputFileName,"read");
    TDirectoryFile *dirResults=(TDirectoryFile*)(fInput->FindObjectAny(Form("HMResults")));
    TList *Results;
    dirResults->GetObject(Form("HMResults"),Results);
    TList* tmpFolder=(TList*)Results->FindObject("Particle0_Particle0");
    TH1F* hSE = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle0");;
    TH1F* hME = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle0");;

    hSE->Scale(1./hSE->Integral(),"width");
    hME->Scale(1./hME->Integral(),"width");

    DLM_Histo<double> dhME_MeV;
    dhME_MeV.SetUp(1);
    dhME_MeV.SetUp(0,hME->GetNbinsX(),hME->GetBinLowEdge(1)*1000.,hME->GetXaxis()->GetBinUpEdge(hME->GetNbinsX())*1000.);
    dhME_MeV.Initialize();

    TH1F* hSE_MeV = new TH1F("hSE_MeV","hSE_MeV",hSE->GetNbinsX(),hSE->GetBinLowEdge(1)*1000.,hSE->GetXaxis()->GetBinUpEdge(hSE->GetNbinsX())*1000.);
    for(unsigned uBin=1; uBin<=hSE->GetNbinsX(); uBin++){
        hSE_MeV->SetBinContent(uBin,hSE->GetBinContent(uBin));
        hSE_MeV->SetBinError(uBin,hSE->GetBinError(uBin));
        dhME_MeV.SetBinContent(uBin,hME->GetBinContent(uBin));
        dhME_MeV.SetBinError(uBin,hME->GetBinError(uBin));
    }

    TH1F* hME_MeV = new TH1F("hME_MeV","hME_MeV",hME->GetNbinsX(),hME->GetBinLowEdge(1)*1000.,hME->GetXaxis()->GetBinUpEdge(hME->GetNbinsX())*1000.);
    for(unsigned uBin=1; uBin<=hME->GetNbinsX(); uBin++){
        hME_MeV->SetBinContent(uBin,hME->GetBinContent(uBin));
        hME_MeV->SetBinError(uBin,hME->GetBinError(uBin));
    }

    TH1F* hCk_ME = (TH1F*)hSE->Clone("hCk_ME");
    hCk_ME->Divide(hME);

    TH1F* hCk_ME_MeV = (TH1F*)hSE_MeV->Clone("hCk_ME_MeV");
    hCk_ME_MeV->Divide(hME);


/// FIT SE WITH A STRETCHED ME

    TF1* fit_Gauss = new TF1("fit_Gauss",GaussSourceTF1,kMin/1000.,kMax/1000.,1);
    fit_Gauss->SetParameter(0,0.5);
    DLM_CleverMcLevyReso CleverMcLevyReso;
    CleverMcLevyReso.InitStability(21,1,2);
    CleverMcLevyReso.InitScale(38,0.1,1.0);
    CleverMcLevyReso.InitRad(513,0,32);
    CleverMcLevyReso.InitType(2);
    CleverMcLevyReso.InitNumMcIter(100000);
    TF1* fit_Levy = new TF1("fit_Levy",&CleverMcLevyReso,&DLM_CleverMcLevyReso::RootEval,kMin/1000.,kMax/1000.,2,"DLM_CleverMcLevyReso","RootEval");
    fit_Levy->SetParameter(0,0.5);
    fit_Levy->SetParLimits(0,0.1,1.0);
    fit_Levy->SetParameter(1,1.5);
    fit_Levy->SetParLimits(1,1.0,2.0);

    fit_Levy->FixParameter(0,0.5);
    fit_Levy->FixParameter(1,1.7);

    hSE->Fit(fit_Gauss,"S, N, R, M");
    TGraph gfSE_Gauss;
    gfSE_Gauss.SetName("gfSE_Gauss");
    gfSE_Gauss.Set(hSE->GetNbinsX());
    for(unsigned uBin=0; uBin<hSE->GetNbinsX(); uBin++){
        double Momentum = hSE->GetBinCenter(uBin+1);
        gfSE_Gauss.SetPoint(uBin,Momentum,fit_Gauss->Eval(Momentum));
    }

    hSE->Fit(fit_Levy,"S, N, R, M");
    TGraph gfSE_Levy;
    gfSE_Levy.SetName("gfSE_Levy");
    gfSE_Levy.Set(hSE->GetNbinsX());
    for(unsigned uBin=0; uBin<hSE->GetNbinsX(); uBin++){
        double Momentum = hSE->GetBinCenter(uBin+1);
        gfSE_Levy.SetPoint(uBin,Momentum,fit_Levy->Eval(Momentum));
    }

    HME_FIT = &dhME_MeV;
    TF1* fitStretchedME = new TF1("fitStretchedME",Fit_StretchME,dhME_MeV.GetLowEdge(0),dhME_MeV.GetUpEdge(0),1);
    //TF1* fitStretchedME = new TF1("fitStretchedME",Fit_StretchME,250,750,1);
    fitStretchedME->SetParameter(0,1);
    fitStretchedME->SetParLimits(0,0.5,1.5);
    printf("hSE->Fit(fitStretchedME,S, N, R, M,,750,2000);\n");
    hSE_MeV->Fit(fitStretchedME,"S, N, R, M","",200,500);
    //hSE_MeV->Fit(fitStretchedME,"S, N, R, M");

    hME->Fit(fit_Gauss,"S, N, R, M");
    TGraph gfME_Gauss;
    gfME_Gauss.SetName("gfME_Gauss");
    gfME_Gauss.Set(hME->GetNbinsX());
    for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
        double Momentum = hME->GetBinCenter(uBin+1);
        gfME_Gauss.SetPoint(uBin,Momentum,fit_Gauss->Eval(Momentum));
    }

    hME->Fit(fit_Levy,"S, N, R, M");
    TGraph gfME_Levy;
    gfME_Levy.SetName("gfME_Levy");
    gfME_Levy.Set(hME->GetNbinsX());
    for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
        double Momentum = hME->GetBinCenter(uBin+1);
        gfME_Levy.SetPoint(uBin,Momentum,fit_Levy->Eval(Momentum));
    }


    double* MomBins_pp = NULL;
    double* FitRegion_pp = NULL;
    unsigned NumMomBins_pp;
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    AnalysisObject.SetUpBinning_pp("pp13TeV_HM_March19",NumMomBins_pp,MomBins_pp,FitRegion_pp,0,0);
    FitRegion_pp[2] = FitRegion_pp[1];
    FitRegion_pp[3] = kMax;
    CATS AB_pp;
    DLM_Ck* Ck_pp;
    AB_pp.SetMomBins(NumMomBins_pp,MomBins_pp);
    AB_pp.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18","Gauss");
    AB_pp.SetAnaSource(0,1.3);
    Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
    Ck_pp->SetSourcePar(0,1.3);
    Ck_pp->SetCutOff(kMax,500);
    Ck_pp->Update();
    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pp,MomBins_pp);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
    AB_pL.SetAnaSource(0,1.3);
    AB_pL.SetNotifications(CATS::nWarning);
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
    Ck_pL->SetSourcePar(0,1.3);
    Ck_pL->Update();
    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_March19","pp");
    //TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_March19","pLambda");
    TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
    DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);
    double lam_pp[5];
    double lam_pL[5];
    AnalysisObject.SetUpLambdaPars_pp("pp13TeV_HM_March19",0,lam_pp);
    AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_March19",0,0,lam_pL);
    CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);
    CkDec_pL.AddContribution(0,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(1,lam_pL[4],DLM_CkDecomposition::cFake);//0.03
    FIT_CK = &CkDec_pp;
    const int NumKnots = 9;
    const double FirstSteps = 6;

    double* Nodes_x = new double [NumKnots];
    Nodes_x[0] = kMin;
    Nodes_x[1] = kMin+FirstSteps;
    Nodes_x[2] = kMin+2*FirstSteps;
    Nodes_x[NumKnots-1] = kMax;
    double NodeLength = (Nodes_x[NumKnots-1]-Nodes_x[1])/double(NumKnots-3);
    printf("NodeLength=%f\n",NodeLength);
    for(int iKnot=3; iKnot<NumKnots-1; iKnot++){
        Nodes_x[iKnot] = Nodes_x[iKnot-1]+NodeLength;
        printf("Nodes_x[%i]=%f\n",iKnot,Nodes_x[iKnot]);
    }

    //Nodes_x[0] = kMin;
    //Nodes_x[NumKnots-1] = kMax;
    //NodeLength = (Nodes_x[NumKnots-1]-Nodes_x[1])/double(NumKnots-1);
    //for(int iKnot=1; iKnot<NumKnots-1; iKnot++){
    //    Nodes_x[iKnot] = Nodes_x[iKnot-1]+NodeLength;
    //    printf("Nodes_x[%i]=%f\n",iKnot,Nodes_x[iKnot]);
    //}




    TF1* fitFemtoSplines_SE = new TF1("fitFemtoSplines_SE",Fit_SE_Splines,kMin,kMax,2+3+NumKnots*2);
    fitFemtoSplines_SE->SetParameter(0,1.3);
    fitFemtoSplines_SE->SetParLimits(0,1.0,1.5);
    fitFemtoSplines_SE->FixParameter(1,2.0);
    fitFemtoSplines_SE->FixParameter(2,NumKnots);
    int KnotCount=0;
    fitFemtoSplines_SE->FixParameter(3+KnotCount++,0);
    //fitFemtoSplines_SE->FixParameter(3+KnotCount++,-4.73e-5);
    //fitFemtoSplines_SE->FixParameter(3+KnotCount++,0);
    fitFemtoSplines_SE->SetParameter(3+KnotCount++,0);
    fitFemtoSplines_SE->SetParLimits(3+KnotCount++,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hSE_MeV->GetBinContent(hSE_MeV->FindBin(Nodes_x[uKnot]));
        fitFemtoSplines_SE->FixParameter(3+2+uKnot,Nodes_x[uKnot]);
        fitFemtoSplines_SE->SetParameter(3+2+NumKnots+uKnot,HistVal);
        fitFemtoSplines_SE->SetParLimits(3+2+NumKnots+uKnot,0,HistVal*2);
        //fitFemtoSplines_SE->SetParameter(3+2+NumKnots+uKnot,0.4);
        //fitFemtoSplines_SE->SetParLimits(3+2+NumKnots+uKnot,0.0,1.2);
    }

    hSE_MeV->Fit(fitFemtoSplines_SE,"S, N, R, M");
    fitFemtoSplines_SE->SetNpx(512);
    double Chi2=fitFemtoSplines_SE->GetChisquare();
    double Ndf=fitFemtoSplines_SE->GetNDF();
    printf("Chi2/ndf=%.2f\n",Chi2/Ndf);
    printf("nSigma=%.2f\n",sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2,int(Ndf))));

    const unsigned FirstBin = hSE_MeV->FindBin(kMin);
    const unsigned LastBin = hSE_MeV->FindBin(kMax-0.01);
    const unsigned NumFitBins = LastBin-FirstBin+1;
    TH1F* hCkExp = new TH1F("hCkExp","hCkExp",NumFitBins,kMin,kMax);
    TGraph gCkTh;
    gCkTh.SetName("gCkTh");
    gCkTh.Set(NumFitBins);
    TGraph gCkTh_ME;
    gCkTh_ME.SetName("gCkTh_ME");
    gCkTh_ME.Set(NumFitBins);
    TGraph gCkTh_ME_Spline;
    gCkTh_ME_Spline.SetName("gCkTh_ME_Spline");
    gCkTh_ME_Spline.Set(NumFitBins);

    TGraph gBl;
    gBl.SetName("gBl");
    gBl.Set(NumFitBins);
    TH1F* hBl_ME = new TH1F("hBl_ME","hBl_ME",NumFitBins,kMin,kMax);
    TH1F* hBl_Res = new TH1F("hBl_Res","hBl_Res",NumFitBins,kMin,kMax);
    TH1F* hBlSpline_Res = new TH1F("hBlSpline_Res","hBlSpline_Res",NumFitBins,kMin,kMax);
    TGraph gBl_Spline;
    gBl_Spline.SetName("gBl_Spline");
    gBl_Spline.Set(NumFitBins);

    for(unsigned uBin=0; uBin<NumFitBins; uBin++){
        double Momentum = hSE_MeV->GetBinCenter(uBin+1);
        hCkExp->SetBinContent(uBin+1,hSE_MeV->GetBinContent(uBin+1)/fitFemtoSplines_SE->Eval(Momentum)*CkDec_pp.EvalCk(Momentum));
        hCkExp->SetBinError(uBin+1,hSE_MeV->GetBinError(uBin+1)/fitFemtoSplines_SE->Eval(Momentum)*CkDec_pp.EvalCk(Momentum));
        gCkTh.SetPoint(uBin,Momentum,CkDec_pp.EvalCk(Momentum));
        gBl.SetPoint(uBin,Momentum,fitFemtoSplines_SE->Eval(Momentum)/CkDec_pp.EvalCk(Momentum));
        hBl_ME->SetBinContent(uBin+1,fitFemtoSplines_SE->Eval(Momentum)/CkDec_pp.EvalCk(Momentum)/hME_MeV->GetBinContent(uBin+1));
        //hBl_ME->SetBinError(uBin+1,fitFemtoSplines_SE->Eval(Momentum)/CkDec_pp.EvalCk(Momentum)/hME_MeV->GetBinContent(uBin+1));
    }


    TF1* fitFemto_Ck = new TF1("fitFemto_Ck",Fit_Ck_Norm,0,300,2+1);
    fitFemto_Ck->SetParameter(0,1.3);
    fitFemto_Ck->SetParLimits(0,1.0,1.5);
    fitFemto_Ck->FixParameter(1,2.0);
    fitFemto_Ck->SetParameter(2,1);
    fitFemto_Ck->SetParLimits(2,0.5,1.5);

    hCk_ME_MeV->Fit(fitFemto_Ck,"S, N, R, M");
    fitFemto_Ck->SetNpx(512);
    double Chi2_Ck=fitFemto_Ck->GetChisquare();
    double Ndf_Ck=fitFemto_Ck->GetNDF();
    printf("Chi2_Ck/ndf_Ck=%.2f\n",Chi2_Ck/Ndf_Ck);
    printf("nSigma_Ck=%.2f\n",sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_Ck,int(Ndf_Ck))));

    for(unsigned uBin=0; uBin<NumFitBins; uBin++){
        double Momentum = hSE_MeV->GetBinCenter(uBin+1);
        gCkTh_ME.SetPoint(uBin,Momentum,CkDec_pp.EvalCk(Momentum));
        hBl_Res->SetBinContent(uBin+1,hCk_ME_MeV->GetBinContent(uBin+1)/fitFemto_Ck->Eval(Momentum)*fitFemto_Ck->GetParameter(2));
        hBl_Res->SetBinError(uBin+1,hCk_ME_MeV->GetBinError(uBin+1)/fitFemto_Ck->Eval(Momentum)*fitFemto_Ck->GetParameter(2));
    }


    FIT_CK = &CkDec_pp;
    const int NumKnots_Ck = 4;
    const double FirstStepsCk = 6;

    double* NodesCk_x = new double [NumKnots];
    NodesCk_x[0] = kMin;
    NodesCk_x[NumKnots_Ck-1] = kMax;
    double NodeLength_Ck = (NodesCk_x[NumKnots-1]-NodesCk_x[1])/double(NumKnots_Ck-1);
    printf("NodeLength_Ck=%f\n",NodeLength_Ck);
    for(int iKnot=1; iKnot<NumKnots_Ck-1; iKnot++){
        NodesCk_x[iKnot] = NodesCk_x[iKnot-1]+NodeLength;
        printf("NodesCk_x[%i]=%f\n",iKnot,NodesCk_x[iKnot]);
    }

    TF1* fitFemtoSpline_Ck = new TF1("fitFemtoSpline_Ck",Fit_SE_Splines,kMin,kMax,2+3+NumKnots*2);
    fitFemtoSpline_Ck->SetParameter(0,1.3);
    fitFemtoSpline_Ck->SetParLimits(0,1.0,1.5);
    fitFemtoSpline_Ck->FixParameter(1,2.0);
    fitFemtoSpline_Ck->FixParameter(2,NumKnots_Ck);
    KnotCount=0;
    fitFemtoSpline_Ck->FixParameter(3+KnotCount++,0);
    //fitFemtoSpline_Ck->FixParameter(3+KnotCount++,-4.73e-5);
    //fitFemtoSpline_Ck->FixParameter(3+KnotCount++,0);
    fitFemtoSpline_Ck->SetParameter(3+KnotCount++,0);
    fitFemtoSpline_Ck->SetParLimits(3+KnotCount++,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots_Ck; uKnot++){
        double HistVal = hCk_ME_MeV->GetBinContent(hCk_ME_MeV->FindBin(NodesCk_x[uKnot]));
        fitFemtoSpline_Ck->FixParameter(3+2+uKnot,NodesCk_x[uKnot]);
        fitFemtoSpline_Ck->SetParameter(3+2+NumKnots_Ck+uKnot,HistVal);
        fitFemtoSpline_Ck->SetParLimits(3+2+NumKnots_Ck+uKnot,0,HistVal*2);
        //fitFemtoSpline_Ck->SetParameter(3+2+NumKnots+uKnot,0.4);
        //fitFemtoSpline_Ck->SetParLimits(3+2+NumKnots+uKnot,0.0,1.2);
    }
    //!make sure the first knot is the same as the 0th one
    fitFemtoSpline_Ck->FixParameter(3+2+NumKnots_Ck+1,1e6);

    hCk_ME_MeV->Fit(fitFemtoSpline_Ck,"S, N, R, M");
    fitFemtoSpline_Ck->SetNpx(512);
    double Chi2_CkSpline=fitFemtoSpline_Ck->GetChisquare();
    double Ndf_CkSpline=fitFemtoSpline_Ck->GetNDF();
    printf("Chi2_CkSpline/Ndf_CkSpline=%.2f\n",Chi2_CkSpline/Ndf_CkSpline);
    printf("nSigma_Ck=%.2f\n",sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_CkSpline,int(Ndf_CkSpline))));

    for(unsigned uBin=0; uBin<NumFitBins; uBin++){
        double Momentum = hSE_MeV->GetBinCenter(uBin+1);
        gCkTh_ME_Spline.SetPoint(uBin,Momentum,CkDec_pp.EvalCk(Momentum));
        hBlSpline_Res->SetBinContent(uBin+1,hCk_ME_MeV->GetBinContent(uBin+1)/fitFemtoSpline_Ck->Eval(Momentum));
        hBlSpline_Res->SetBinError(uBin+1,hCk_ME_MeV->GetBinError(uBin+1)/fitFemtoSpline_Ck->Eval(Momentum));
        gBl_Spline.SetPoint(uBin,Momentum,fitFemtoSpline_Ck->Eval(Momentum)/CkDec_pp.EvalCk(Momentum));
        //hBl_Spline->SetBinError(uBin+1,hCk_ME_MeV->GetBinError(uBin+1)/fitFemto_Ck->Eval(Momentum)*fitFemto_Ck->GetParameter(2));
    }


//TF1("myfunc",myfunction,0,10,2);
    TFile* fOutput = new TFile(OutputFolderName+"fOutput.root","recreate");
    hSE_MeV->Write();
    fitFemtoSplines_SE->Write();
    hME_MeV->Write();
    hCk_ME_MeV->Write();
    fitFemto_Ck->Write();
    fitFemtoSpline_Ck->Write();
    hBl_ME->Write();
    hBl_Res->Write();
    hBlSpline_Res->Write();
    gBl.Write();
    gBl_Spline.Write();
    hCkExp->Write();
    gCkTh.Write();
    gCkTh_ME.Write();
    gCkTh_ME_Spline.Write();
    hSE->Write();
    gfSE_Gauss.Write();
    gfSE_Levy.Write();
    hME->Write();
    gfME_Gauss.Write();
    gfME_Levy.Write();
    fitStretchedME->Write();
    hCk_ME->Write();

    delete hSE_MeV;
    delete hME_MeV;
    delete hCk_ME_MeV;
    delete hCk_ME;
    delete hCkExp;
    delete hBl_ME;
    delete hBl_Res;
    delete fit_Gauss;
    delete fit_Levy;
    delete fitStretchedME;
    delete fitFemtoSplines_SE;
    delete fitFemto_Ck;
    delete fInput;
    delete fOutput;
    delete [] MomBins_pp;
    delete [] FitRegion_pp;
}

void Fit_pL_Splines(){
    //The 1D Boltzmann is actually Gauss with var = s^2 = kT/m
    //=> the distribution of relative momenta between particles following a Boltzmann
    //have the same PDF as S(r), only that we have k as an argument! Isn't that Great!
    //=> I should be able to fit the N(k) with the function we use for the Gaussian source!
    // and if it does not fit, why not use a Levy :D

    const TString FittingMode_pL = "Spline3_5";
    const TString InputFileName = "/home/dmihaylov/CernBox/HM13TeV/AnalysisData/ClosePairRej/SelectedPairs/AnalysisResults.root";
    const TString OutputFolderName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Fit_pL_Splines/";;

    const double kMin = 0;
    const double kMax = 540;

    TFile* fInput = new TFile(InputFileName,"read");
    TDirectoryFile *dirResults=(TDirectoryFile*)(fInput->FindObjectAny(Form("HMResults")));
    TList *Results;
    dirResults->GetObject(Form("HMResults"),Results);
    TList* tmpFolder=(TList*)Results->FindObject("Particle0_Particle2");
    TH1F* hSE = (TH1F*)tmpFolder->FindObject("SEDist_Particle0_Particle2");;
    TH1F* hME = (TH1F*)tmpFolder->FindObject("MEDist_Particle0_Particle2");;

    hSE->Rebin(3);
    hME->Rebin(3);
    hSE->Scale(1./hSE->Integral()/3.,"width");
    hME->Scale(1./hME->Integral()/3.,"width");

    TH1F* hSE_MeV = new TH1F("hSE_MeV","hSE_MeV",hSE->GetNbinsX(),hSE->GetBinLowEdge(1)*1000.,hSE->GetXaxis()->GetBinUpEdge(hSE->GetNbinsX())*1000.);
    for(unsigned uBin=1; uBin<=hSE->GetNbinsX(); uBin++){
        hSE_MeV->SetBinContent(uBin,hSE->GetBinContent(uBin));
        hSE_MeV->SetBinError(uBin,hSE->GetBinError(uBin));
    }

    TH1F* hME_MeV = new TH1F("hME_MeV","hME_MeV",hME->GetNbinsX(),hME->GetBinLowEdge(1)*1000.,hME->GetXaxis()->GetBinUpEdge(hME->GetNbinsX())*1000.);
    for(unsigned uBin=1; uBin<=hME->GetNbinsX(); uBin++){
        hME_MeV->SetBinContent(uBin,hME->GetBinContent(uBin));
        hME_MeV->SetBinError(uBin,hME->GetBinError(uBin));
    }

    TH1F* hCk_ME = (TH1F*)hSE->Clone("hCk_ME");
    hCk_ME->Divide(hME);


    TF1* fit_Gauss = new TF1("fit_Gauss",GaussSourceTF1,kMin/1000.,kMax/1000.,1);
    fit_Gauss->SetParameter(0,0.5);
    DLM_CleverMcLevyReso CleverMcLevyReso;
    CleverMcLevyReso.InitStability(21,1,2);
    CleverMcLevyReso.InitScale(38,0.1,1.0);
    CleverMcLevyReso.InitRad(513,0,32);
    CleverMcLevyReso.InitType(2);
    CleverMcLevyReso.InitNumMcIter(100000);
    TF1* fit_Levy = new TF1("fit_Levy",&CleverMcLevyReso,&DLM_CleverMcLevyReso::RootEval,kMin/1000.,kMax/1000.,2,"DLM_CleverMcLevyReso","RootEval");
    fit_Levy->SetParameter(0,0.5);
    fit_Levy->SetParLimits(0,0.1,1.0);
    fit_Levy->SetParameter(1,1.5);
    fit_Levy->SetParLimits(1,1.0,2.0);

    fit_Levy->FixParameter(0,0.5);
    fit_Levy->FixParameter(1,1.7);

    hSE->Fit(fit_Gauss,"S, N, R, M");
    TGraph gfSE_Gauss;
    gfSE_Gauss.SetName("gfSE_Gauss");
    gfSE_Gauss.Set(hSE->GetNbinsX());
    for(unsigned uBin=0; uBin<hSE->GetNbinsX(); uBin++){
        double Momentum = hSE->GetBinCenter(uBin+1);
        gfSE_Gauss.SetPoint(uBin,Momentum,fit_Gauss->Eval(Momentum));
    }

    hSE->Fit(fit_Levy,"S, N, R, M");
    TGraph gfSE_Levy;
    gfSE_Levy.SetName("gfSE_Levy");
    gfSE_Levy.Set(hSE->GetNbinsX());
    for(unsigned uBin=0; uBin<hSE->GetNbinsX(); uBin++){
        double Momentum = hSE->GetBinCenter(uBin+1);
        gfSE_Levy.SetPoint(uBin,Momentum,fit_Levy->Eval(Momentum));
    }

    hME->Fit(fit_Gauss,"S, N, R, M");
    TGraph gfME_Gauss;
    gfME_Gauss.SetName("gfME_Gauss");
    gfME_Gauss.Set(hME->GetNbinsX());
    for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
        double Momentum = hME->GetBinCenter(uBin+1);
        gfME_Gauss.SetPoint(uBin,Momentum,fit_Gauss->Eval(Momentum));
    }

    hME->Fit(fit_Levy,"S, N, R, M");
    TGraph gfME_Levy;
    gfME_Levy.SetName("gfME_Levy");
    gfME_Levy.Set(hME->GetNbinsX());
    for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
        double Momentum = hME->GetBinCenter(uBin+1);
        gfME_Levy.SetPoint(uBin,Momentum,fit_Levy->Eval(Momentum));
    }


    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    unsigned NumMomBins_pL;
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    AnalysisObject.SetUpBinning_pL("pp13TeV_HM_March19",NumMomBins_pL,MomBins_pL,FitRegion_pL,0,0);
    FitRegion_pL[2] = FitRegion_pL[1];
    FitRegion_pL[3] = kMax;
    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
    AB_pL.SetAnaSource(0,1.3);
    AB_pL.SetNotifications(CATS::nWarning);
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
    Ck_pL->SetSourcePar(0,1.3);
    Ck_pL->Update();
    //TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_March19","pp");
    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_March19","pLambda");
    //TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
    //DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);
    //double lam_pp[5];
    double lam_pL[5];
    //AnalysisObject.SetUpLambdaPars_pp("pp13TeV_HM_March19",0,lam_pp);
    AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_March19",0,0,lam_pL);
    //CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    //CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
    //CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);
    CkDec_pL.AddContribution(0,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(1,lam_pL[4],DLM_CkDecomposition::cFake);//0.03
    FIT_CK = &CkDec_pL;
    const int NumKnots = 7;
    const double FirstSteps = 18;

    double* Nodes_x = new double [NumKnots];
    Nodes_x[0] = kMin;
    Nodes_x[1] = kMin+FirstSteps;
    Nodes_x[NumKnots-1] = kMax;
    double NodeLength = (Nodes_x[NumKnots-1]-Nodes_x[1])/double(NumKnots-2);
    printf("NodeLength=%f\n",NodeLength);
    for(int iKnot=2; iKnot<NumKnots-1; iKnot++){
        Nodes_x[iKnot] = Nodes_x[iKnot-1]+NodeLength;
        printf("Nodes_x[%i]=%f\n",iKnot,Nodes_x[iKnot]);
    }


    Nodes_x[0] = 0;
    Nodes_x[1] = 30;
    Nodes_x[2] = 140;
    Nodes_x[3] = 250;
    Nodes_x[4] = 360;
    Nodes_x[5] = 450;
    Nodes_x[6] = 540;
    //Nodes_x[7] = 540;
    //Nodes_x[8] = 540;
    //Nodes_x[8] = 540;

    TF1* fitFemtoSplines_SE = new TF1("fitFemtoSplines_SE",Fit_SE_Splines,kMin,kMax,2+3+NumKnots*2);
    fitFemtoSplines_SE->SetParameter(0,1.3);
    fitFemtoSplines_SE->SetParLimits(0,1.0,1.8);
    fitFemtoSplines_SE->FixParameter(1,2.0);
    fitFemtoSplines_SE->FixParameter(2,NumKnots);
    int KnotCount=0;
    fitFemtoSplines_SE->FixParameter(3+KnotCount++,0);
    //fitFemtoSplines_SE->FixParameter(3+KnotCount++,-4.73e-5);
    //fitFemtoSplines_SE->FixParameter(3+KnotCount++,0);
    fitFemtoSplines_SE->SetParameter(3+KnotCount++,0);
    fitFemtoSplines_SE->SetParLimits(3+KnotCount++,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hSE_MeV->GetBinContent(hSE_MeV->FindBin(Nodes_x[uKnot]));
        fitFemtoSplines_SE->FixParameter(3+2+uKnot,Nodes_x[uKnot]);
        fitFemtoSplines_SE->SetParameter(3+2+NumKnots+uKnot,HistVal);
        fitFemtoSplines_SE->SetParLimits(3+2+NumKnots+uKnot,0,HistVal*2);
        //fitFemtoSplines_SE->SetParameter(3+2+NumKnots+uKnot,0.4);
        //fitFemtoSplines_SE->SetParLimits(3+2+NumKnots+uKnot,0.0,1.2);
    }

    hSE_MeV->Fit(fitFemtoSplines_SE,"S, N, R, M");
    fitFemtoSplines_SE->SetNpx(512);
    double Chi2=fitFemtoSplines_SE->GetChisquare();
    double Ndf=fitFemtoSplines_SE->GetNDF();
    printf("Chi2/ndf=%.2f\n",Chi2/Ndf);
    printf("nSigma=%.2f\n",sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2,int(Ndf))));


    const unsigned FirstBin = hSE_MeV->FindBin(kMin);
    const unsigned LastBin = hSE_MeV->FindBin(kMax-0.01);
    const unsigned NumFitBins = LastBin-FirstBin+1;
    TH1F* hCkExp = new TH1F("hCkExp","hCkExp",NumFitBins,kMin,kMax);
    TGraph gCkTh;
    gCkTh.SetName("gCkTh");
    gCkTh.Set(NumFitBins);
    TGraph gBl;
    gBl.SetName("gBl");
    gBl.Set(NumFitBins);
    TH1F* hBl_ME = new TH1F("hBl_ME","hBl_ME",NumFitBins,kMin,kMax);
    for(unsigned uBin=0; uBin<NumFitBins; uBin++){
        double Momentum = hSE_MeV->GetBinCenter(uBin+1);
        hCkExp->SetBinContent(uBin+1,hSE_MeV->GetBinContent(uBin+1)/fitFemtoSplines_SE->Eval(Momentum)*CkDec_pL.EvalCk(Momentum));
        hCkExp->SetBinError(uBin+1,hSE_MeV->GetBinError(uBin+1)/fitFemtoSplines_SE->Eval(Momentum)*CkDec_pL.EvalCk(Momentum));
        gCkTh.SetPoint(uBin,Momentum,CkDec_pL.EvalCk(Momentum));
        gBl.SetPoint(uBin,Momentum,fitFemtoSplines_SE->Eval(Momentum)/CkDec_pL.EvalCk(Momentum));
        hBl_ME->SetBinContent(uBin+1,fitFemtoSplines_SE->Eval(Momentum)/CkDec_pL.EvalCk(Momentum)/hME_MeV->GetBinContent(uBin+1));
        //hBl_ME->SetBinError(uBin+1,fitFemtoSplines_SE->Eval(Momentum)/CkDec_pp.EvalCk(Momentum)/hME_MeV->GetBinContent(uBin+1));
    }

//TF1("myfunc",myfunction,0,10,2);
    TFile* fOutput = new TFile(OutputFolderName+"fOutput.root","recreate");
    hSE_MeV->Write();
    hME_MeV->Write();
    hBl_ME->Write();
    fitFemtoSplines_SE->Write();
    gBl.Write();
    hCkExp->Write();
    gCkTh.Write();
    hSE->Write();
    gfSE_Gauss.Write();
    gfSE_Levy.Write();
    hME->Write();
    gfME_Gauss.Write();
    gfME_Levy.Write();
    hCk_ME->Write();

    delete hSE_MeV;
    delete hME_MeV;
    delete hCk_ME;
    delete hCkExp;
    delete hBl_ME;
    delete fit_Gauss;
    delete fit_Levy;
    delete fitFemtoSplines_SE;
    delete fInput;
    delete fOutput;
    delete [] MomBins_pL;
    delete [] FitRegion_pL;
}


void FitTheSameEventOnly(){

}

//as a way to see what the bl might be
void RatioBetweenLevy(){

    const double Size1 = 360;
    const double Stab1 = 1.54;

    const double Size2 = 400;
    const double Stab2 = 1.46;

    double kMin = 0;
    double kMax = 4000;

    DLM_CleverMcLevyReso CleverMcLevyReso;
    CleverMcLevyReso.InitStability(21,1.45,1.55);
    CleverMcLevyReso.InitScale(38,300,500);
    CleverMcLevyReso.InitRad(257,0,4000);
    CleverMcLevyReso.InitType(2);
    CleverMcLevyReso.InitNumMcIter(20000000);
    TF1* fit_Levy = new TF1("fit_Levy",&CleverMcLevyReso,&DLM_CleverMcLevyReso::RootEval,kMin,kMax,2,"DLM_CleverMcLevyReso","RootEval");
    fit_Levy->FixParameter(0,Size1);
    fit_Levy->FixParameter(1,Stab1);

    TGraph Lev1;
    Lev1.SetName("Lev1");
    TGraph Lev2;
    Lev2.SetName("Lev2");
    TGraph LevRatio;
    LevRatio.SetName("LevRatio");

    const unsigned NumBins = 256;
    const double BinLen = (kMax-kMin)/double(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        double Momentum = kMin+BinLen*0.5+BinLen*double(uBin);
        Lev1.SetPoint(uBin,Momentum,fit_Levy->Eval(Momentum));
    }

    fit_Levy->FixParameter(0,Size2);
    fit_Levy->FixParameter(1,Stab2);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        double Momentum = kMin+BinLen*0.5+BinLen*double(uBin);
        double OldVal, Dummy;
        Lev1.GetPoint(uBin,Dummy,OldVal);
        double NewVal = fit_Levy->Eval(Momentum);
        Lev2.SetPoint(uBin,Momentum,NewVal);
        if(OldVal>1e-6) LevRatio.SetPoint(uBin,Momentum,NewVal/OldVal);
        else LevRatio.SetPoint(uBin,Momentum,0);
    }

    TFile* fOut = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/RatioBetweenLevy/fOut.root","recreate");
    Lev1.Write();
    Lev2.Write();
    LevRatio.Write();

    delete fOut;

}


TH1F* CDS_HADD(const unsigned& uDist, const CatsEvent* KittyEvent, TH1F* histo){
  return NULL;
    /*
    TH1F* hAdd = NULL;
    const int NumMomBins = histo->GetNbinsX();
    const double kMin = histo->GetBinLowEdge(1);
    const double kMax = histo->GetYaxis()->GetBinUpEdge(NumMomBins);

    const double EtaMax = 0.8;
    const double pT_Min = 0.4;
    const double pT_Max = 4.0;
    const unsigned NumPart1 = KittyEvent->GetNumParticles1();
    const unsigned NumPart2 = KittyEvent->GetNumParticles2();
    if(!NumPart1||!NumPart2) return NULL;
    CatsEvent CurrentEvent(KittyEvent->GetParticleType1(0)->GetPid(),KittyEvent->GetParticleType2(0)->GetPid());
    const unsigned NumPairs = KittyEvent->GetNumPairs();
    const unsigned NumPart1 = KittyEvent->GetNumParticles1();
    const unsigned NumPart2 = KittyEvent->GetNumParticles2();
    double pT_1,pT_2,eta_1,eta_2;
    if(uDist<5){
        hAdd = new TH1F("hAdd","hAdd",NumMomBins,kMin,kMax);
    }
    else{
        return NULL;
    }
    if(uDist==0){
        for(unsigned uPair=0; uPair<NumPairs; uPair++){
            CatsParticle& Part1 = *KittyEvent->GetParticlePair(0)->GetParticle(0);
            CatsParticle& Part2 = *KittyEvent->GetParticlePair(0)->GetParticle(1);
            if(Part1.GetPt()<pT_Min||Part1.GetPt()>pT_Max) continue;
            if(Part2.GetPt()<pT_Min||Part2.GetPt()>pT_Max) continue;
            hAdd->Fill(KittyEvent->GetParticlePair(0)->GetP()*0.5);
        }
    }
    else if(uDist==1){
            CatsParticle& Part1 = *KittyEvent->GetParticlePair(0)->GetParticle(0);
            CatsParticle& Part2 = *KittyEvent->GetParticlePair(0)->GetParticle(1);
            for(unsigned uPart2=0; uPart2<NumPart2; uPart2++){

            }


            CatsParticle& Part1 = *KittyEvent->GetParticlePair(0)->GetParticle(0);
            CatsParticle& Part2 = *KittyEvent->GetParticlePair(0)->GetParticle(1);
            CurrentEvent.AddParticle();

            Part1 = KittyEvent->GetParticlePair(0)->GetParticle(0);
            Part2 = KittyEvent->GetParticlePair(0)->GetParticle(1);
            if(Part1->GetPt()<pT_Min||Part1->GetPt()>pT_Max) continue;
            if(Part2->GetPt()<pT_Min||Part2->GetPt()>pT_Max) continue;
            hAdd->Fill(KittyEvent->GetParticlePair(0)->GetP()*0.5);
    }
    else if(uDist==2){

    }
    else if(uDist==3){

    }
    else if(uDist==4){

    }

    hAdd->Sumw2();
    histo->Add(hAdd);

    delete hAdd;
    */
}

void CompareReferenceSamples(const TString DataSetDescr){
    const double kMin = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMax = DataSetDescr=="pp"?4:DataSetDescr=="pLambda"?4:DataSetDescr=="pXim"?4:4;
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/CompareReferenceSamples/";
    const unsigned NumMomBins = DataSetDescr=="pp"?200:DataSetDescr=="pLambda"?100:DataSetDescr=="pXim"?100:100;
    const TString OutFileBaseName = OutputFolder+DataSetDescr;
    const TString InputFileName = DataSetDescr=="pp"?TransportFile_pp_Alice:DataSetDescr=="pLambda"?TransportFile_pL_Alice:
        DataSetDescr=="pXim"?TString("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_100KiLrz_pXim.f19"):"";
    const unsigned NumBlankHeaderLines=3;
    const unsigned MaxPairsToRead = 4e6;
    const unsigned HighMultLimit = 128;
    const unsigned MixingDepth = 8;
    const unsigned RandomDepth = 16;
    const double DimiSmearPhi = (90.)*DegToRad;
    int pdgID[2] = {2212, DataSetDescr=="pp"?2212:DataSetDescr=="pLambda"?3122:DataSetDescr=="pXim"?3312:0};

    const unsigned NumDistr = 5;
    TString* DistrName = new TString [NumDistr];
    TString* DistrDescr = new TString [NumDistr];
    DistrName[0] = "EPOS Signal";
    DistrDescr[0] = "SIG";
    DistrName[1] = "EPOS SE";
    DistrDescr[1] = "SE";
    DistrName[2] = "Random Phi";
    DistrDescr[2] = "RndPhi";
    DistrName[3] = "Dimi's method";
    DistrDescr[3] = "DimiPhi";
    DistrName[4] = "Stravinsky method";
    DistrDescr[4] = "Strav";

    TH1F** hSE = new TH1F* [NumDistr];
    for(unsigned uDist=0; uDist<NumDistr; uDist++){
        hSE[uDist] = new TH1F("h"+DistrDescr[uDist],"h"+DistrDescr[uDist],NumMomBins,kMin,kMax);
    }

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    unsigned NumTotalPairs=0;
    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;
    unsigned NumSePairs;
    unsigned NumMePairs;


    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;

    while(!feof(InFile)){
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
//printf("TotNumEvents=%u\n",TotNumEvents);
        TotNumEvents++;
        if(NumPartInEvent>int(HighMultLimit)) RejectedHighMultEvents++;
        CatsParticle KittyParticle;
        CatsEvent* KittyEvent = new CatsEvent(pdgID[0],pdgID[1]);
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            if(NumTotalPairs>=MaxPairsToRead) continue;
            if(NumPartInEvent>int(HighMultLimit)) continue;

//ALICE ACCEPTANCE
if(KittyParticle.GetP()<0.4) continue;
if(fabs(KittyParticle.GetPseudoRap())>0.8) continue;//!

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            KittyEvent->AddParticle(KittyParticle);
        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)

        KittyEvent->ComputeParticlePairs();

        NumTotalPairs += KittyEvent->GetNumPairs();

        for(unsigned uDist=0; uDist<NumDistr; uDist++){
            CDS_HADD(uDist,KittyEvent,hSE[uDist]);
            //if(hADD){hSE[uDist]->Add(hADD); delete hADD; hADD=NULL;}
        }

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }
        delete KittyEvent;
    }//while(!feof(InFile))

    if(ProgressBar){
        printf("\r\033[K");
    }

    delete [] DistrDescr;
    delete [] DistrName;
    for(unsigned uDist=0; uDist<NumDistr; uDist++){
        delete hSE[uDist];
    }
    delete [] hSE;
}

void SE_for_Vale(const TString DataSetDescr){
    const double kMin = 0;
    const double kMax = 4;
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/SE_for_Vale/";
    const unsigned NumMomBins = 200;
    const TString OutFileBaseName = OutputFolder+DataSetDescr;
    const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19";
    const unsigned NumBlankHeaderLines = 3;
    const unsigned MaxPairsToRead = 1e9;
    const unsigned HighMultLimit = 128;
    int pdgID[2] = {2212,3122};

    TH1F* hSE = new TH1F("hSE","hSE",NumMomBins,kMin,kMax);

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    unsigned NumTotalPairs=0;
    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;

    while(!feof(InFile)){
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
        TotNumEvents++;
        if(NumPartInEvent>int(HighMultLimit)) RejectedHighMultEvents++;
        CatsParticle KittyParticle;
        CatsEvent* KittyEvent = new CatsEvent(pdgID[0],pdgID[1]);
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            if(NumTotalPairs>=MaxPairsToRead) continue;
            if(NumPartInEvent>int(HighMultLimit)) continue;

            ///ALICE ACCEPTANCE
            //if(KittyParticle.GetP()<0.4) continue;
            if(fabs(KittyParticle.GetPseudoRap())>0.8) continue;

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            KittyEvent->AddParticle(KittyParticle);
        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)

        KittyEvent->ComputeParticlePairs();

        NumTotalPairs += KittyEvent->GetNumPairs();

        for(unsigned uPair=0; uPair<KittyEvent->GetNumPairs(); uPair++){
            hSE->Fill((KittyEvent->GetParticlePair(uPair).GetP())*0.5);
        }

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }
        delete KittyEvent;
    }//while(!feof(InFile))

    if(ProgressBar){
        printf("\r\033[K");
    }

    hSE->Sumw2();

    TFile OutputFile(OutFileBaseName+"_SE.root","recreate");
    hSE->Write();

    delete hSE;
    delete cdummy;
}


//the original version, where we assume that the coordinate systems (CS) is fixed from one primordial protons
//there is no separation between primordial-reso and reso-reso case
void ReferenceSampleStudy_1(const TString& TranModDescr, const TString& DataSetDescr){
    const double kMin = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMax = DataSetDescr=="pp"?4000:DataSetDescr=="pLambda"?4000:DataSetDescr=="pXim"?4000:4000;
    const double kMinZoom = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMaxZoom = DataSetDescr=="pp"?0.6:DataSetDescr=="pLambda"?0.6:DataSetDescr=="pXim"?0.6:0.6;
    //const double kNormMin = 0.400;
    //const double kNormMax = 0.600;
    //const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/ReferenceSampleStudy_1/";
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/";
    const unsigned NumMomBins = DataSetDescr=="pp"?400:DataSetDescr=="pLambda"?200:DataSetDescr=="pXim"?200:200;
    const TString OutFileBaseName = OutputFolder+TranModDescr+"_"+DataSetDescr;
    //const TString InputFileName = DataSetDescr=="pp"?TransportFile_pp_Alice:DataSetDescr=="pLambda"?TransportFile_pL_Alice:
    //    DataSetDescr=="pXim"?TString("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_100KiLrz_pXim.f19"):"";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/ThermalSpectrum10_pp2/thermalTest2_PRS=1.30_PTS=1.30.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_prim_030219.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pXim_PRIM_4PI.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pOmega_PRIM_4PI.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pReso_PRIM_4PI.f19";
    const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp80/pp80_pReso_PRIM_4PI.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp80/pp80_LamReso_PRIM_4PI.f19";

    const unsigned NumBlankHeaderLines=3;
    const unsigned MaxPairsToRead = 200e6;
    const unsigned HighMultLimit = 128;
    const double SmearPhi = (45.)*DegToRad;
    int pdgID[2] = {DataSetDescr=="LambdaLambda"?3122:2212,
    DataSetDescr=="pp"?2212:DataSetDescr=="pLambda"?3122:DataSetDescr=="LambdaLambda"?3122:DataSetDescr=="pXim"?3312:DataSetDescr=="pOmega"?3334:0};

    if(DataSetDescr=="pReso"){
        pdgID[0] = 1120;
        pdgID[1] = 0;
    }
    else if(DataSetDescr=="LamReso"){
        pdgID[0] = 2130;
        pdgID[1] = 0;
    }
    else if(DataSetDescr=="pResoReso"||DataSetDescr=="LamResoReso"){
        pdgID[0] = 0;
        pdgID[1] = 0;
    }


//pdgID[0]=3122;
//pdgID[1]=2212;
    unsigned NumTotalPairs=0;
    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;
    //0 is the same event
    //1 is mixed event without multiplicity weighting
    //2 is mixed event with multiplicity weighting
    //3 is rotated phi
    //4 is is rotated phi with conservation
    const unsigned NumMethods = 5;
    TString MethodShortDescr[NumMethods];
    MethodShortDescr[0] = "SE";
    MethodShortDescr[1] = "MEraw";
    MethodShortDescr[2] = "MEmult";
    MethodShortDescr[3] = "RotPhi";
    MethodShortDescr[4] = "RotPhiCons";
    unsigned MixingDepth[NumMethods];
    MixingDepth[0] = 1;
    MixingDepth[1] = 8;
    MixingDepth[2] = 8;
    MixingDepth[3] = 8;
    MixingDepth[4] = 8;

    TH1F** hDist = new TH1F* [NumMethods];
    TH1F** hCk = new TH1F* [NumMethods];
    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        hDist[uMeth] = new TH1F(TString::Format("h%s",MethodShortDescr[uMeth].Data()),TString::Format("h%s",MethodShortDescr[uMeth].Data()),NumMomBins,kMin,kMax);
        hCk[uMeth] = new TH1F(TString::Format("hCk%s",MethodShortDescr[uMeth].Data()),TString::Format("hCk%s",MethodShortDescr[uMeth].Data()),NumMomBins,kMin,kMax);
    }

    //the pair is boosted to the rest frame of the 1st particle.
    //the angle between r* and k* is plotted vs the momentum of the 2nd particle
    //this is needed in order to investigate the assumption of the resonance source, where r* and k* are taken as parallel
    //As no EPOS files with resonances were available, this is best applied to the case of a primordial source, since all primordial
    //particles should have similar kinematics.
    ////if we are in the double reso case (both PID==0), than we decay one of the resonances into a proton+pi or Lam+pi and use the daughter to
    ////fix the coordinate system. Than we look at angle between the two particles in THAT system
    TH2F* h_rkAngle_Mom2 = new TH2F("h_rkAngle_Mom2","h_rkAngle_Mom2",64,0,Pi,64,0,1024);
    //the angle between the momenta of the two primordial particles
    TH2F* h_kkAngle_Mom2 = new TH2F("h_kkAngle_Mom2","h_kkAngle_Mom2",64,0,Pi,64,0,1024);

    CatsParticle KittyParticle;
    CatsEvent*** KittyEvent = new CatsEvent** [NumMethods];
    CatsEvent*** KittyFilteredEvent = new CatsEvent** [NumMethods];
    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        KittyEvent[uMeth] = new CatsEvent* [MixingDepth[uMeth]];
        KittyFilteredEvent[uMeth] = new CatsEvent* [MixingDepth[uMeth]];
        for(unsigned uMix=0; uMix<MixingDepth[uMeth]; uMix++){
            KittyEvent[uMeth][uMix] = new CatsEvent(pdgID[0],pdgID[1]);
            KittyFilteredEvent[uMeth][uMix] = new CatsEvent(pdgID[0],pdgID[1]);
        }
    }

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;
    TRandom3 rangen(11);
    while(!feof(InFile)){
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
        TotNumEvents++;
//printf("TotNumEvents=%u\n",TotNumEvents);
        if(NumPartInEvent>int(HighMultLimit)) RejectedHighMultEvents++;
//printf("NumPartInEvent=%u\n",NumPartInEvent);
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            int ParticlePID = KittyParticle.GetPid();
            if(DataSetDescr.Contains("Reso")&&KittyParticle.GetPid()!=pdgID[0])KittyParticle.SetPid(0);

            switch(ParticlePID){
                //PDG
                case 32224 : KittyParticle.SetWidth(0.35); break;//Delta 1600
                case 32214 : KittyParticle.SetWidth(0.35); break;
                case 32114 : KittyParticle.SetWidth(0.35); break;
                case 2226 : KittyParticle.SetWidth(0.33); break;//Delta 1905
                case 2126 : KittyParticle.SetWidth(0.33); break;
                case 1216 : KittyParticle.SetWidth(0.33); break;
                case 12112 : KittyParticle.SetWidth(0.30); break;//N 1440
                case 12212 : KittyParticle.SetWidth(0.30); break;
                case 12224 : KittyParticle.SetWidth(0.30); break;//Delta 1700
                case 12214 : KittyParticle.SetWidth(0.30); break;
                case 12114 : KittyParticle.SetWidth(0.30); break;
                case 22222 : KittyParticle.SetWidth(0.25); break;//Delta 1910
                case 22122 : KittyParticle.SetWidth(0.25); break;
                case 21212 : KittyParticle.SetWidth(0.25); break;
                case 31214 : KittyParticle.SetWidth(0.20); break;//N 1440
                case 32124 : KittyParticle.SetWidth(0.20); break;
                case 32112 : KittyParticle.SetWidth(0.165); break;//N 1650
                case 32212 : KittyParticle.SetWidth(0.165); break;
                case 22212 : KittyParticle.SetWidth(0.15); break;//N 1535
                case 22112 : KittyParticle.SetWidth(0.15); break;
                case 2116 : KittyParticle.SetWidth(0.15); break;//N 1675
                case 2216 : KittyParticle.SetWidth(0.15); break;
                case 2222 : KittyParticle.SetWidth(0.145); break;//Delta 1620
                case 2122 : KittyParticle.SetWidth(0.145); break;
                case 1212 : KittyParticle.SetWidth(0.145); break;
                case 12116 : KittyParticle.SetWidth(0.13); break;//N 1680
                case 12216 : KittyParticle.SetWidth(0.13); break;
                case 2224 : KittyParticle.SetWidth(0.118); break;//Delta 1230
                case 2214 : KittyParticle.SetWidth(0.118); break;
                case 2114 : KittyParticle.SetWidth(0.118); break;
                case 1214 : KittyParticle.SetWidth(0.115); break;//N 1520
                case 2124 : KittyParticle.SetWidth(0.115); break;
                case 42212 : KittyParticle.SetWidth(0.10); break;//N 1710
                case 42112 : KittyParticle.SetWidth(0.10); break;
                case 21214 : KittyParticle.SetWidth(0.10); break;//N 1700
                case 22124 : KittyParticle.SetWidth(0.10); break;
                //EPOS
                case 1111 : KittyParticle.SetWidth(0.118); break;//Delta 1230
                case 1121 : KittyParticle.SetWidth(0.118); break;
                case 1221 : KittyParticle.SetWidth(0.118); break;
                case 1122 : KittyParticle.SetWidth(0.20); break;//N 1440
                case 1222 : KittyParticle.SetWidth(0.20); break;
                case 1123 : KittyParticle.SetWidth(0.15); break;//N 1530 -> 1535
                case 1223 : KittyParticle.SetWidth(0.15); break;
                case 1112 : KittyParticle.SetWidth(0.145); break;//Delta 1620
                case 1124 : KittyParticle.SetWidth(0.145); break;
                case 1224 : KittyParticle.SetWidth(0.145); break;
                case 1125 : KittyParticle.SetWidth(0.15); break;//N 1665 -> 1675
                case 1225 : KittyParticle.SetWidth(0.15); break;
                case 1113 : KittyParticle.SetWidth(0.30); break;//Delta 1700
                case 1126 : KittyParticle.SetWidth(0.30); break;
                case 1226 : KittyParticle.SetWidth(0.30); break;
                case 1127 : KittyParticle.SetWidth(0.10); break;//N 1710
                case 1227 : KittyParticle.SetWidth(0.10); break;
                case 1114 : KittyParticle.SetWidth(0.25); break;//Delta 1925 -> 1910
                case 1128 : KittyParticle.SetWidth(0.25); break;
                case 1228 : KittyParticle.SetWidth(0.25); break;
                default : KittyParticle.SetWidth(0); break;
            }


            //continue;
            //ALICE ACCEPTANCE
            //if(KittyParticle.GetP()<0.4) continue;
            //if(fabs(KittyParticle.GetPseudoRap())>0.8) continue;//!

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
                for(unsigned uMix=0; uMix<MixingDepth[uMeth]; uMix++){
                    KittyEvent[uMeth][uMix]->AddParticle(KittyParticle);
                }
            }
        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)
        for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
            for(unsigned uMix=0; uMix<MixingDepth[uMeth]; uMix++){
                if(uMeth==0){
                    //do nothing
                }
                else if(uMeth==1){
                    //mixed events
                }
                else if(uMeth==2){
                    //mixed events with mult.
                }
                else if(uMeth==3){
                    //rotated phi
                    KittyEvent[uMeth][uMix]->RandomizeMomentumPhi(SmearPhi,"","Uniform");
                }
                else if(uMeth==4){
                    //rotated phi with E conservation
                    KittyEvent[uMeth][uMix]->RandomizeMomentumPhi(SmearPhi,"EnergyConservation","Uniform");
                }
                for(unsigned uPart=0; uPart<KittyEvent[uMeth][uMix]->GetNumParticles1(); uPart++){
                    //if(KittyEvent[uMeth][uMix]->GetParticleType1(uPart).GetP()<0.4) continue;
                    //if(fabs(KittyEvent[uMeth][uMix]->GetParticleType1(uPart).GetPseudoRap())>0.8) continue;
                    KittyFilteredEvent[uMeth][uMix]->AddParticle(KittyEvent[uMeth][uMix]->GetParticleType1(uPart));
                }
                for(unsigned uPart=0; uPart<KittyEvent[uMeth][uMix]->GetNumParticles2(); uPart++){
                    //avoid double counting for identical particles
                    if(KittyEvent[uMeth][uMix]->GetSameType()) break;
                    //if(KittyEvent[uMeth][uMix]->GetParticleType2(uPart).GetP()<0.4) continue;
                    //if(fabs(KittyEvent[uMeth][uMix]->GetParticleType2(uPart).GetPseudoRap())>0.8) continue;
                    KittyFilteredEvent[uMeth][uMix]->AddParticle(KittyEvent[uMeth][uMix]->GetParticleType2(uPart));
                }
                KittyFilteredEvent[uMeth][uMix]->ComputeParticlePairs();
            }
        }

        NumTotalPairs+=KittyFilteredEvent[0][0]->GetNumPairs();

        for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
            for(unsigned uMix=0; uMix<MixingDepth[uMeth]; uMix++){
                for(unsigned uPair=0; uPair<KittyFilteredEvent[uMeth][uMix]->GetNumPairs(); uPair++){
                    int pid0 = KittyFilteredEvent[uMeth][uMix]->GetParticlePair(uPair).GetParticle(0).GetPid();
                    int pid1 = KittyFilteredEvent[uMeth][uMix]->GetParticlePair(uPair).GetParticle(1).GetPid();
                    //make sure we have particles of the desired PID
                    if((pdgID[0]!=pid0||pdgID[1]!=pid1) && (pdgID[0]!=pid1||pdgID[1]!=pid0)) continue;

                    hDist[uMeth]->Fill(KittyFilteredEvent[uMeth][uMix]->GetParticlePair(uPair).GetP()*1000.);
                    //fill the angular distributions only for the default same event sample
                    if(uMeth==0&&uMix==0){
                        CatsParticle Particle1 = KittyFilteredEvent[uMeth][uMix]->GetParticlePair(uPair).GetParticle(0);
                        CatsParticle Particle2 = KittyFilteredEvent[uMeth][uMix]->GetParticlePair(uPair).GetParticle(1);

                        TVector3 vP1;
                        vP1.SetXYZ(Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
                        TVector3 vP2;
                        vP2.SetXYZ(Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
                        double Angle_P1P2 = vP1.Angle(vP2);
                        h_kkAngle_Mom2->Fill(Angle_P1P2,(vP1-vP2).Mag()*1000.);

                        //printf("------------------\n");
                        //printf("p1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
                        //printf("p2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
                        //boost into the rest frame of the first particle
                        CatsLorentzVector BoostVector;
                        BoostVector.Set(0,0,0,0,Particle1.GetE(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
                        //printf("BV = (%.3f,%.3f,%.3f)\n",BoostVector.GetPx(),BoostVector.GetPy(),BoostVector.GetPz());
                        Particle1.Boost(BoostVector);
                        Particle2.Boost(BoostVector);
                        TVector3 vecR;
                        vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
                        TVector3 vecK;
                        vecK.SetXYZ(Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
                        double Angle = vecR.Angle(vecK);
                        if(Angle<0) Angle += Pi;
                        if(Angle>Pi) Angle -= Pi;

                        //printf("p1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
                        //printf("p2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
                        //printf("Angle = %f\n",Angle*RadToDeg);

                        h_rkAngle_Mom2->Fill(Angle,Particle2.GetP()*1000.);

                    }
                }
            }

        }

        for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
            for(unsigned uMix=0; uMix<MixingDepth[uMeth]; uMix++){
                KittyEvent[uMeth][uMix]->Reset();
                KittyFilteredEvent[uMeth][uMix]->Reset();
            }
        }

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }

    }//while(!feof(InFile))

    if(ProgressBar){
        printf("\r\033[K");
    }

    //normalize
    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        hDist[uMeth]->Sumw2();
        hDist[uMeth]->Scale(1./hDist[uMeth]->Integral(),"width");
    }

    //compute Ck
    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            hCk[uMeth]->SetBinContent(uBin,hDist[0]->GetBinContent(uBin));
            hCk[uMeth]->SetBinError(uBin,hDist[0]->GetBinError(uBin));
        }
        hCk[uMeth]->Divide(hDist[uMeth]);
    }

    TFile* fOut = new TFile(OutFileBaseName+".root","recreate");
    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        hDist[uMeth]->Write();
    }
    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        hCk[uMeth]->Write();
    }

    h_rkAngle_Mom2->GetXaxis()->SetTitle("Angle (rad)");
    h_rkAngle_Mom2->GetYaxis()->SetTitle("p_{res} (MeV)");

    h_rkAngle_Mom2->Write();
    h_kkAngle_Mom2->Write();

    for(unsigned uMeth=0; uMeth<NumMethods; uMeth++){
        for(unsigned uMix=0; uMix<MixingDepth[uMeth]; uMix++){
            delete KittyEvent[uMeth][uMix];
        }
        delete [] KittyEvent[uMeth];
        delete hDist[uMeth];
    }
    delete [] KittyEvent;
    delete [] hDist;
    delete h_rkAngle_Mom2;
    delete h_kkAngle_Mom2;
    delete fOut;
}



void ResoInfo(const int& ResoPid, int& ParentPid, double& Width, const TString& System){
  //pion resonances
  if(System=="pi_Sig"){
    if(ResoPid==211||ResoPid==120){Width=0; ParentPid=0;}
    else if(ResoPid==3122||ResoPid==2130){Width=0; ParentPid=0;}
    else {Width=0; ParentPid=1;}
    return;
  }
  if(System=="piReso_Sig"||System=="pi_SigReso"||System=="pi_Sig"||System=="piReso_SigReso"){
    //primordial
    if(abs(ResoPid)==120||abs(ResoPid)==211||abs(ResoPid)==2130||abs(ResoPid)==3122) {Width=0; ParentPid=0;}
    else{
      switch (abs(ResoPid)) {
        //pion resonances
        case 120: Width=0; ParentPid=0; break;
        case 111: Width=1; ParentPid=211; break;
        case 121: Width=1; ParentPid=211; break;
        case 131: Width=1; ParentPid=211; break;
        case 231: Width=1; ParentPid=211; break;
        case 10111: Width=1; ParentPid=211; break;
        case 1111: Width=1; ParentPid=211; break;
        case 2221: Width=1; ParentPid=211; break;
        case 1221: Width=1; ParentPid=211; break;
        case 1121: Width=1; ParentPid=211; break;
        case 1122: Width=1; ParentPid=211; break;
        case 1125: Width=1; ParentPid=211; break;
        case 1126: Width=1; ParentPid=211; break;
        case 1113: Width=1; ParentPid=211; break;
        case 1127: Width=1; ParentPid=211; break;
        case 1128: Width=1; ParentPid=211; break;
        case 1114: Width=1; ParentPid=211; break;

        case 211: Width=0; ParentPid=0; break;
        case 113: Width=1; ParentPid=211; break;
        case 213: Width=1; ParentPid=211; break;
        case 323: Width=1; ParentPid=211; break;
        case 313: Width=1; ParentPid=211; break;
        case 20113: Width=1; ParentPid=211; break;
        case 2224: Width=1; ParentPid=211; break;
        //case 1114: Width=1; ParentPid=211; break;//duplicates, but in both cases it is a pion reso, so okay
        case 2114: Width=1; ParentPid=211; break;
        case 2214: Width=1; ParentPid=211; break;
        case 12112: Width=1; ParentPid=211; break;
        case 22212: Width=1; ParentPid=211; break;
        case 2122: Width=1; ParentPid=211; break;
        case 2222: Width=1; ParentPid=211; break;
        case 32212: Width=1; ParentPid=211; break;
        case 12214: Width=1; ParentPid=211; break;
        case 12224: Width=1; ParentPid=211; break;

        //sigma resonances -> well actually its Lambda used as a dummy!
        case 2130: Width=0; ParentPid=0; break;
        case 1131: Width=1; ParentPid=3122; break;
        case 1231: Width=1; ParentPid=3122; break;
        case 2231: Width=1; ParentPid=3122; break;
        case 1132: Width=1; ParentPid=3122; break;
        case 1236: Width=1; ParentPid=3122; break;
        case 2232: Width=1; ParentPid=3122; break;
        case 1133: Width=1; ParentPid=3122; break;
        case 1237: Width=1; ParentPid=3122; break;
        case 2233: Width=1; ParentPid=3122; break;
        case 1134: Width=1; ParentPid=3122; break;
        case 1239: Width=1; ParentPid=3122; break;
        case 2234: Width=1; ParentPid=3122; break;

        case 3122: Width=0; ParentPid=0; break;
        case 3224: Width=1; ParentPid=3122; break;
        case 3214: Width=1; ParentPid=3122; break;
        case 3114: Width=1; ParentPid=3122; break;
        case 13124: Width=1; ParentPid=3122; break;
        case 13222: Width=1; ParentPid=3122; break;
        case 13212: Width=1; ParentPid=3122; break;
        case 13112: Width=1; ParentPid=3122; break;
        case 13224: Width=1; ParentPid=3122; break;
        case 13214: Width=1; ParentPid=3122; break;
        case 13114: Width=1; ParentPid=3122; break;
        case 3216: Width=1; ParentPid=3122; break;
        case 3226: Width=1; ParentPid=3122; break;
        case 3116: Width=1; ParentPid=3122; break;
        case 33122: Width=1; ParentPid=3122; break;
        case 13226: Width=1; ParentPid=3122; break;
        case 13216: Width=1; ParentPid=3122; break;
        case 13116: Width=1; ParentPid=3122; break;
        case 13324: Width=1; ParentPid=3122; break;
        case 13314: Width=1; ParentPid=3122; break;
        case 3126: Width=1; ParentPid=3122; break;
        case 23222: Width=1; ParentPid=3122; break;
        case 23212: Width=1; ParentPid=3122; break;
        case 23112: Width=1; ParentPid=3122; break;
        default : printf("\033[1;33mWARNING:\033[0m Unexpected pi_Sig resonances!\n"); Width=0; ParentPid=1; break;
      }
    }
  }
  //
  else if(System=="p_KaonReso"||System=="pReso_KaonReso"||System=="pReso_Kaon"){
    switch (abs(ResoPid)) {
      case 1120: Width=0; ParentPid=0; break;
      case 1111: Width=1; ParentPid=2212; break;
      case 1121: Width=1; ParentPid=2212; break;
      case 1221: Width=1; ParentPid=2212; break;
      case 1122: Width=1; ParentPid=2212; break;
      case 1222: Width=1; ParentPid=2212; break;
      case 1123: Width=1; ParentPid=2212; break;
      case 1223: Width=1; ParentPid=2212; break;
      case 1112: Width=1; ParentPid=2212; break;
      case 1124: Width=1; ParentPid=2212; break;
      case 1224: Width=1; ParentPid=2212; break;
      case 1125: Width=1; ParentPid=2212; break;
      case 1225: Width=1; ParentPid=2212; break;
      case 1113: Width=1; ParentPid=2212; break;
      case 1126: Width=1; ParentPid=2212; break;
      case 1226: Width=1; ParentPid=2212; break;
      case 1127: Width=1; ParentPid=2212; break;
      case 1227: Width=1; ParentPid=2212; break;
      case 1114: Width=1; ParentPid=2212; break;
      case 1128: Width=1; ParentPid=2212; break;
      case 1228: Width=1; ParentPid=2212; break;

      case 2212: Width=0; ParentPid=0; break;
      case 2224: Width=1; ParentPid=2212; break;
      case 2214: Width=1; ParentPid=2212; break;
      case 2114: Width=1; ParentPid=2212; break;
      case 12112: Width=1; ParentPid=2212; break;
      case 1214: Width=1; ParentPid=2212; break;
      case 12116: Width=1; ParentPid=2212; break;
      case 22212: Width=1; ParentPid=2212; break;
      case 12212: Width=1; ParentPid=2212; break;
      case 2124: Width=1; ParentPid=2212; break;
      case 2226: Width=1; ParentPid=2212; break;
      case 12224: Width=1; ParentPid=2212; break;
      case 2116: Width=1; ParentPid=2212; break;
      case 31214: Width=1; ParentPid=2212; break;
      case 12216: Width=1; ParentPid=2212; break;
      case 12214: Width=1; ParentPid=2212; break;
      case 22112: Width=1; ParentPid=2212; break;
      case 2126: Width=1; ParentPid=2212; break;
      case 32224: Width=1; ParentPid=2212; break;
      case 32112: Width=1; ParentPid=2212; break;
      case 2216: Width=1; ParentPid=2212; break;
      case 2222: Width=1; ParentPid=2212; break;
      case 32124: Width=1; ParentPid=2212; break;
      case 32214: Width=1; ParentPid=2212; break;
      case 32212: Width=1; ParentPid=2212; break;
      case 12114: Width=1; ParentPid=2212; break;
      case 2122: Width=1; ParentPid=2212; break;
      case 1216: Width=1; ParentPid=2212; break;
      case 22222: Width=1; ParentPid=2212; break;
      case 32114: Width=1; ParentPid=2212; break;
      case 42212: Width=1; ParentPid=2212; break;
      case 22122: Width=1; ParentPid=2212; break;
      case 1212: Width=1; ParentPid=2212; break;
      case 42112: Width=1; ParentPid=2212; break;
      case 21214: Width=1; ParentPid=2212; break;
      case 21212: Width=1; ParentPid=2212; break;
      case 22124: Width=1; ParentPid=2212; break;

      case 130: Width=0; ParentPid=0; break;
      case 231: Width=1; ParentPid=321; break;
      case 131: Width=1; ParentPid=321; break;
      case 1220: Width=1; ParentPid=321; break;
      case 2130: Width=1; ParentPid=321; break;
      case 1233: Width=1; ParentPid=321; break;
      case 1234: Width=1; ParentPid=321; break;
      case 2330: Width=1; ParentPid=321; break;
      case 1330: Width=1; ParentPid=321; break;
      case 3331: Width=1; ParentPid=321; break;
      case 2234: Width=1; ParentPid=321; break;
      case 1134: Width=1; ParentPid=321; break;
      case 1239: Width=1; ParentPid=321; break;

      case 321: Width=0; ParentPid=0; break;
      case 313: Width=1; ParentPid=321; break;
      case 323: Width=1; ParentPid=321; break;
      case 2112: Width=1; ParentPid=321; break;
      case 315: Width=1; ParentPid=321; break;
      case 10313: Width=1; ParentPid=321; break;
      case 10323: Width=1; ParentPid=321; break;
      case 9000225: Width=1; ParentPid=321; break;
      case 325: Width=1; ParentPid=321; break;
      case 9010221: Width=1; ParentPid=321; break;
      case 9000111: Width=1; ParentPid=321; break;
      case 10311: Width=1; ParentPid=321; break;
      case 335: Width=1; ParentPid=321; break;
      case 10333: Width=1; ParentPid=321; break;
      case 3124: Width=1; ParentPid=321; break;
      case 215: Width=1; ParentPid=321; break;
      case 9000223: Width=1; ParentPid=321; break;
      case 10321: Width=1; ParentPid=321; break;
      case 30313: Width=1; ParentPid=321; break;
      case 20333: Width=1; ParentPid=321; break;
      case 327: Width=1; ParentPid=321; break;
      case 100215: Width=1; ParentPid=321; break;
      case 20223: Width=1; ParentPid=321; break;
      case 225: Width=1; ParentPid=321; break;
      default: printf("\033[1;33mWARNING:\033[0m Unexpected p_Kaon resonances!\n"); Width=0; ParentPid=1; break;
    }
  }
  else{
    switch(ResoPid){
        //PDG
        case 32224 : Width = 0.35; ParentPid = 2212; break;//Delta 1600
        case 32214 : Width = 0.35; ParentPid = 2212; break;
        case 32114 : Width = 0.35; ParentPid = 2212; break;
        case 2226 : Width = 0.33; ParentPid = 2212; break;//Delta 1905
        case 2126 : Width = 0.33; ParentPid = 2212; break;
        case 1216 : Width = 0.33; ParentPid = 2212; break;
        case 12112 : Width = 0.30; ParentPid = 2212; break;//N 1440
        case 12212 : Width = 0.30; ParentPid = 2212; break;
        case 12224 : Width = 0.30; ParentPid = 2212; break;//Delta 1700
        case 12214 : Width = 0.30; ParentPid = 2212; break;
        case 12114 : Width = 0.30; ParentPid = 2212; break;
        case 22222 : Width = 0.25; ParentPid = 2212; break;//Delta 1910
        case 22122 : Width = 0.25; ParentPid = 2212; break;
        case 21212 : Width = 0.25; ParentPid = 2212; break;
        case 31214 : Width = 0.20; ParentPid = 2212; break;//N 1440
        case 32124 : Width = 0.20; ParentPid = 2212; break;
        case 32112 : Width = 0.165; ParentPid = 2212; break;//N 1650
        case 32212 : Width = 0.165; ParentPid = 2212; break;
        case 22212 : Width = 0.15; ParentPid = 2212; break;//N 1535
        case 22112 : Width = 0.15; ParentPid = 2212; break;
        case 2116 : Width = 0.15; ParentPid = 2212; break;//N 1675
        case 2216 : Width = 0.15; ParentPid = 2212; break;
        case 2222 : Width = 0.145; ParentPid = 2212; break;//Delta 1620
        case 2122 : Width = 0.145; ParentPid = 2212; break;
        case 1212 : Width = 0.145; ParentPid = 2212; break;
        case 12116 : Width = 0.13; ParentPid = 2212; break;//N 1680
        case 12216 : Width = 0.13; ParentPid = 2212; break;
        case 2224 : Width = 0.118; ParentPid = 2212; break;//Delta 1230
        case 2214 : Width = 0.118; ParentPid = 2212; break;
        case 2114 : Width = 0.118; ParentPid = 2212; break;
        case 1214 : Width = 0.115; ParentPid = 2212; break;//N 1520
        case 2124 : Width = 0.115; ParentPid = 2212; break;
        case 42212 : Width = 0.10; ParentPid = 2212; break;//N 1710
        case 42112 : Width = 0.10; ParentPid = 2212; break;
        case 21214 : Width = 0.10; ParentPid = 2212; break;//N 1700
        case 22124 : Width = 0.10; ParentPid = 2212; break;
        //for lambda
        case 3224 : Width = 0.0358; ParentPid = 3122; break;//S*+
        case 3214 : Width = 0.036; ParentPid = 3122; break;//S*0
        case 3114 : Width = 0.0394; ParentPid = 3122; break;//S*-
        case 13124 : Width = 0.06; ParentPid = 3122; break;//Lambda 1690
        case 13222 : Width = 0.1; ParentPid = 3122; break;//S+(1660)
        case 13212 : Width = 0.1; ParentPid = 3122; break;//S0(1660)
        case 13112 : Width = 0.1; ParentPid = 3122; break;//S-(1660)
        case 13224 : Width = 0.06; ParentPid = 3122; break;//S+(1670)
        case 13214 : Width = 0.06; ParentPid = 3122; break;//S0(1670)
        case 13114 : Width = 0.06; ParentPid = 3122; break;//S-(1670)
        case 3216 : Width = 0.12; ParentPid = 3122; break;//S+(1775)
        case 3226 : Width = 0.12; ParentPid = 3122; break;//S0(1775)
        case 3116 : Width = 0.12; ParentPid = 3122; break;//S-(1775)
        case 33122 : Width = 0.035; ParentPid = 3122; break;//Lambda 1670
        case 13226 : Width = 0.12; ParentPid = 3122; break;//S+(1915)
        case 13216 : Width = 0.12; ParentPid = 3122; break;//S0(1915)
        case 13116 : Width = 0.12; ParentPid = 3122; break;//S-(1915)
        case 13324 : Width = 0.024; ParentPid = 3122; break;//Xi0 1820
        case 13314 : Width = 0.024; ParentPid = 3122; break;//Xi- 1820
        case 3126 : Width = 0.08; ParentPid = 3122; break;//Lambda 1820
        case 23222 : Width = 0.09; ParentPid = 3122; break;//S+(1750)
        case 23212 : Width = 0.09; ParentPid = 3122; break;//S0(1750)
        case 23112 : Width = 0.09; ParentPid = 3122; break;//S-(1750)
        //case 42212 : Width = 0.1; break;//N+ 1710
        //case 42112 : Width = 0.1; break;//N0 1710
        //case 32124 : Width = 0.2; break;//N+ 1720
        //case 31214 : Width = 0.2; break;//N0 1720
        //case 32212 : Width = 0.165; break;//N+ 1650
        //case 32112 : Width = 0.165; break;//N0 1650

        //EPOS
        case 1111 : Width = 0.118; ParentPid = 1120; break;//Delta 1230
        case 1121 : Width = 0.118; ParentPid = 1120; break;
        case 1221 : Width = 0.118; ParentPid = 1120; break;
        case 1122 : Width = 0.20; ParentPid = 1120; break;//N 1440
        case 1222 : Width = 0.20; ParentPid = 1120; break;
        case 1123 : Width = 0.15; ParentPid = 1120; break;//N 1530 -> 1535
        case 1223 : Width = 0.15; ParentPid = 1120; break;
        case 1112 : Width = 0.145; ParentPid = 1120; break;//Delta 1620
        case 1124 : Width = 0.145; ParentPid = 1120; break;
        case 1224 : Width = 0.145; ParentPid = 1120; break;
        case 1125 : Width = 0.15; ParentPid = 1120; break;//N 1665 -> 1675
        case 1225 : Width = 0.15; ParentPid = 1120; break;
        case 1113 : Width = 0.30; ParentPid = 1120; break;//Delta 1700
        case 1126 : Width = 0.30; ParentPid = 1120; break;
        case 1226 : Width = 0.30; ParentPid = 1120; break;
        case 1127 : Width = 0.10; ParentPid = 1120; break;//N 1710
        case 1227 : Width = 0.10; ParentPid = 1120; break;
        case 1114 : Width = 0.25; ParentPid = 1120; break;//Delta 1925 -> 1910
        case 1128 : Width = 0.25; ParentPid = 1120; break;
        case 1228 : Width = 0.25; ParentPid = 1120; break;
        //for lambda
        case 1131 : Width = 0.036; ParentPid = 2130; break;//S1385
        case 1231 : Width = 0.036; ParentPid = 2130; break;
        case 2231 : Width = 0.039; ParentPid = 2130; break;
        case 1132 : Width = 0.08; ParentPid = 2130; break;//S1665 (avg of 1660 and 1670)
        case 1236 : Width = 0.08; ParentPid = 2130; break;
        case 2232 : Width = 0.08; ParentPid = 2130; break;
        case 1133 : Width = 0.12; ParentPid = 2130; break;//S1775
        case 1237 : Width = 0.12; ParentPid = 2130; break;
        case 2233 : Width = 0.12; ParentPid = 2130; break;
        case 1134 : Width = 0.12; ParentPid = 2130; break;
        case 1239 : Width = 0.12; ParentPid = 2130; break;
        case 2234 : Width = 0.12; ParentPid = 2130; break;
        default : Width = 0; ParentPid=0; break;
    }
  }
}


//New version of ReferenceSampleStudy_1, where we now decay all resonances and take the daughters to fix the coordinate system
void ReferenceSampleStudy_2(const TString& TranModDescr, const TString& DataSetDescr){
    //if true, it decays the parents multiple times, and the solution which results in the smallest k* of the daughters is considered
    //the idea: enhance the statistics for the region we care about.
    const bool PickBestDecay = false;
    const double kMin = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMax = DataSetDescr=="pp"?4000:DataSetDescr=="pLambda"?4000:DataSetDescr=="pXim"?4000:4000;
    const double kMinZoom = DataSetDescr=="pp"?0:DataSetDescr=="pLambda"?0:DataSetDescr=="pXim"?0:0;
    const double kMaxZoom = DataSetDescr=="pp"?0.6:DataSetDescr=="pLambda"?0.6:DataSetDescr=="pXim"?0.6:0.6;
    //const double kNormMin = 0.400;
    //const double kNormMax = 0.600;
    //const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/ReferenceSampleStudy_1/";
    //const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_3/";
    //const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/";
    const TString OutputFolder = TString::Format("%s/CatsFiles/Source/EposAngularDist/",GetCernBoxDimi());
    const unsigned NumMomBins = DataSetDescr=="pp"?400:DataSetDescr=="pLambda"?200:DataSetDescr=="pXim"?200:200;
    const TString OutFileBaseName = OutputFolder+TranModDescr+"_"+DataSetDescr;
    //const TString InputFileName = DataSetDescr=="pp"?TransportFile_pp_Alice:DataSetDescr=="pLambda"?TransportFile_pL_Alice:
    //    DataSetDescr=="pXim"?TString("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_100KiLrz_pXim.f19"):"";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/ThermalSpectrum10_pp2/thermalTest2_PRS=1.30_PTS=1.30.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_prim_030219.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pXim_PRIM_4PI.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pOmega_PRIM_4PI.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pReso_PRIM_4PI_VER4.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pReso_4PI_VER5.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pReso_Oct2019_4PI.f19";

    TString InputFileName;
    if( (DataSetDescr.Contains("Xim")||DataSetDescr.Contains("Omega"))&&!DataSetDescr.Contains("Lam"))
        InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_pResoXimOmega_Oct2019_4PI_ReducedWeights.f19";
    else if(DataSetDescr=="Lam_Lam"||DataSetDescr=="Lam_LamReso"||DataSetDescr=="LamReso_LamReso")
        InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp200/pp200_LamResoLamReso_Oct2019_4PI_ReducedWeights.f19";
    else if(DataSetDescr=="pi_pi"||DataSetDescr=="pi_piReso"||DataSetDescr=="piReso_piReso")
        //InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Max/220320.f19";
        //InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Max/270320_with_omega.f19";
        InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Max/300320_width50_omegaYes.root";
        //InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Max/300320_width50_omegaNo.root";
    else if(DataSetDescr=="pi_Sig"||DataSetDescr=="piReso_Sig"||DataSetDescr=="pi_SigReso"||DataSetDescr=="piReso_SigReso"){
      InputFileName = TString::Format("%s/CatsFiles/Source/EposRawOutput/pi_Sig_Dummy.f19",GetCernBoxDimi());
    }
    else if(DataSetDescr=="p_Kaon"||DataSetDescr=="pReso_Kaon"||DataSetDescr=="p_KaonReso"||DataSetDescr=="pReso_KaonReso")
        InputFileName = TString::Format("%s/scratch6/dmihaylov/OutputEPOS/13TeV/EPOS_20200121/pKch_2021_80.f19",GetNx2Folder());
    else if(DataSetDescr=="p_Phi"||DataSetDescr=="pReso_Phi")
        InputFileName = "/home/dimihayl/Mount/nx3/scratch6/dmihaylov/OutputEPOS/13TeV/EPOS_20200121/pPhi_All.f19";
    else if(DataSetDescr=="Lam_Xim"||DataSetDescr=="LamReso_Xim")
      InputFileName = "/home/dimihayl/Mount/nx3/scratch6/dmihaylov/OutputEPOS/13TeV/EPOS_20200121/LamXi_All.f19";
    else
        //InputFileName = "$PATH_NX1/scratch6/dmihaylov/OutputEPOS/13TeV/EPOS_LBF_pp200/pp200_pResoLamReso_Oct2019_4PI_ReducedWeights.f19";
        InputFileName = "/home/dimihayl/Mount/nx1/scratch6/dmihaylov/OutputEPOS/13TeV/EPOS_LBF_pp200/pp200_pResoLamReso_Oct2019_4PI_ReducedWeights.f19";
//
    //AvgResoMass = 1.354190 for VER2, around 1.38 for the first version. Basically VER2 has inly ior==0 reso
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp80/pp80_pReso_PRIM_4PI_VER2.f19";
    //const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp80/pp80_LamReso_PRIM_4PI.f19";

    const unsigned NumBlankHeaderLines=3;
    const unsigned MaxPairsToRead = 2000e6;
    const unsigned HighMultLimit = 128;
    const double SmearPhi = (45.)*DegToRad;
    int pdgID[2] = {DataSetDescr=="LambdaLambda"?3122:2212,
    DataSetDescr=="pp"?2212:DataSetDescr=="pLambda"?3122:DataSetDescr=="LambdaLambda"?3122:DataSetDescr=="pXim"?3312:DataSetDescr=="pOmega"?3334:0};
    int Original_eposID[2];
    int Original_pdgID[2];
    //masses of the decay products (daughters)
    //a maximum of four body decay is allowed. It is assumed that [0] Daughter is the only one of interest
    double DaughterMassP1[4];
    double DaughterMassP2[4];
    double ResoMass[2];

    //[0] = 2 body //[1] = 3 body //[2] = 4 body
    double FractionsNbody[3];
    double FractionsNbodyC[3];//cumulative

    //see the comment at ResoInfo
    bool UseResoInfo = true;
    //the minumum mass of the resonances required to perform a multi body decay
    //the fractions refer only to resonances of above that mass!!!
    double MultiBodyThreshold = 0;

    //p-p correlation
    if(DataSetDescr=="p_p"){
        pdgID[0] = 2212;
        pdgID[1] = 2212;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 1120;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 2212;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="p_pReso"||DataSetDescr=="p_pReso_3body"){
        bool Use_3body = DataSetDescr.Contains("3body");

        pdgID[0] = 2212;
        pdgID[1] = 0;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 1120;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 2212;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_p*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = Use_3body?Mass_pic*0.001:-1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = 1.36;

        //refers only to guys above MultiBodyThreshold
        //these numbers come from the idea of having a total of 1% of all resonances
        //have a 3-body decay, by trial and mistake I found that 2.3% of the resonances
        //with mass above 1400 MeV consition do the job.
        FractionsNbody[0] = Use_3body?0.955:1;
        FractionsNbody[1] = Use_3body?0.045:0;
        FractionsNbody[2] = 0;
        MultiBodyThreshold = 1.4;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pReso_pReso"||DataSetDescr=="pReso_pReso_3body"){
        bool Use_3body = DataSetDescr.Contains("3body");

        pdgID[0] = 0;
        pdgID[1] = 0;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 1120;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 2212;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = Use_3body?Mass_pic*0.001:-1;;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_p*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = Use_3body?Mass_pic*0.001:-1;;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = 1.36;

        //refers only to guys above MultiBodyThreshold
        //these numbers come from the idea of having a total of 1% of all resonances
        //have a 3-body decay, by trial and mistake I found that 2.3% of the resonances
        //with mass above 1400 MeV consition do the job.
        FractionsNbody[0] = Use_3body?0.977:1;
        FractionsNbody[1] = Use_3body?0.023:0;
        FractionsNbody[2] = 0;
        MultiBodyThreshold = 1.4;

        UseResoInfo = true;
    }
    //p-Lambda correlation
    else if(DataSetDescr=="p_LamReso"){
        pdgID[0] = 2212;
        pdgID[1] = 1;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_L*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = 1.46;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
//! WRONG PIDs?
    else if(DataSetDescr=="pReso_Lam"){
        pdgID[0] = 0;
        pdgID[1] = 2130;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = -1;

        pdgID[0] = 0;
        pdgID[1] = 3122;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;//was wrong set to 1.46
        ResoMass[1] = -1;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pReso_LamReso"){
        pdgID[0] = 0;
        pdgID[1] = 1;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_L*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = 1.46;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pReso_Xim"){
        pdgID[0] = 0;
        pdgID[1] = 2330;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = -1;

        pdgID[0] = 0;
        pdgID[1] = 3312;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 2330;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 3312;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pReso_Omega"){
        pdgID[0] = 0;
        pdgID[1] = 3331;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = -1;

        pdgID[0] = 0;
        pdgID[1] = 3334;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 3331;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 3334;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    //LL correlation
    else if(DataSetDescr=="Lam_Lam"){
        pdgID[0] = 3122;
        pdgID[1] = 3122;
        Original_eposID[0] = 2130;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 3122;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="Lam_LamReso"){
        pdgID[0] = 3122;
        pdgID[1] = 0;
        Original_eposID[0] = 2130;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 3122;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_L*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = 1.46;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="LamReso_LamReso"){
        pdgID[0] = 0;
        pdgID[1] = 0;
        Original_eposID[0] = 2130;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 3122;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = Mass_L*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_L*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.46;
        ResoMass[1] = 1.46;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    //pipi correlation
    else if(DataSetDescr=="pi_pi"){
        pdgID[0] = 211;
        pdgID[1] = 211;
        Original_eposID[0] = 120;
        Original_eposID[1] = 120;
        Original_pdgID[0] = 211;
        Original_pdgID[1] = 211;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = -1;

        FractionsNbody[0] = 0.934;
        FractionsNbody[1] = 0.037;
        FractionsNbody[2] = 0.029;

        UseResoInfo = false;
    }
    else if(DataSetDescr=="pi_piReso"){
        pdgID[0] = 211;
        pdgID[1] = 0;
        Original_eposID[0] = 120;
        Original_eposID[1] = 120;
        Original_pdgID[0] = 211;
        Original_pdgID[1] = 211;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_pic*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = Mass_pic*0.001;
        DaughterMassP2[3] = Mass_pic*0.001;

        ResoMass[0] = -1;
        ResoMass[1] = 1.124;

        FractionsNbody[0] = 0.934;
        FractionsNbody[1] = 0.037;
        FractionsNbody[2] = 0.029;

        UseResoInfo = false;
    }
    else if(DataSetDescr=="piReso_piReso"){
        pdgID[0] = 0;
        pdgID[1] = 0;
        Original_eposID[0] = 120;
        Original_eposID[1] = 120;
        Original_pdgID[0] = 211;
        Original_pdgID[1] = 211;

        DaughterMassP1[0] = Mass_pic*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = Mass_pic*0.001;
        DaughterMassP1[3] = Mass_pic*0.001;

        DaughterMassP2[0] = Mass_pic*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = Mass_pic*0.001;
        DaughterMassP2[3] = Mass_pic*0.001;

        ResoMass[0] = 1.124;
        ResoMass[1] = 1.124;

        FractionsNbody[0] = 0.934;
        FractionsNbody[1] = 0.037;
        FractionsNbody[2] = 0.029;

        UseResoInfo = false;
    }
    //Modification to Dimis code
    //piSig correlation for now use Lambda=Sigma, neglect many body-decays
    else if(DataSetDescr=="pi_Sig"){
        pdgID[0] = 211;
        pdgID[1] = 3122;
        Original_eposID[0] = 120;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 211;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="piReso_Sig"){
        pdgID[0] = 0;
        pdgID[1] = 3122;
        Original_eposID[0] = 120;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 211;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = Mass_pic*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = Mass_pic*0.001;
        DaughterMassP1[3] = Mass_pic*0.001;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.124;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pi_SigReso"){
        pdgID[0] = 211;
        pdgID[1] = 1;
        Original_eposID[0] = 120;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 211;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_L*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = 1.46;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="piReso_SigReso"){
        pdgID[0] = 0;
        pdgID[1] = 1;
        Original_eposID[0] = 120;
        Original_eposID[1] = 2130;
        Original_pdgID[0] = 211;
        Original_pdgID[1] = 3122;

        DaughterMassP1[0] = Mass_pic*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = Mass_pic*0.001;
        DaughterMassP1[3] = Mass_pic*0.001;

        DaughterMassP2[0] = Mass_L*0.001;
        DaughterMassP2[1] = Mass_pic*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.124;
        ResoMass[1] = 1.46;

        FractionsNbody[0] = 1;
        FractionsNbody[1] = 0;
        FractionsNbody[2] = 0;

        UseResoInfo = true;
    }

    //pK correlation
    else if(DataSetDescr=="p_Kaon"){
        pdgID[0] = 2212;
        pdgID[1] = 321;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 130;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 321;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="p_KaonReso"){
        pdgID[0] = 2212;
        pdgID[1] = 1;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 130;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 321;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_Kch*0.001;
        DaughterMassP2[1] = Mass_pi0*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = 1.054;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pReso_Kaon"){
        pdgID[0] = 0;
        pdgID[1] = 321;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 130;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 321;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pReso_KaonReso"){
        pdgID[0] = 0;
        pdgID[1] = 1;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 130;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 321;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = Mass_Kch*0.001;
        DaughterMassP2[1] = Mass_pi0*0.001;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = 1.054;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }
    //pPhi correlation
    else if(DataSetDescr=="p_Phi"){
        pdgID[0] = 2212;
        pdgID[1] = 333;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 331;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 333;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="pReso_Phi"){
        pdgID[0] = 0;
        pdgID[1] = 333;
        Original_eposID[0] = 1120;
        Original_eposID[1] = 331;
        Original_pdgID[0] = 2212;
        Original_pdgID[1] = 333;

        DaughterMassP1[0] = Mass_p*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.36;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }
    //LamXim correlation
    else if(DataSetDescr=="Lam_Xim"){
        pdgID[0] = 3122;
        pdgID[1] = 3312;
        Original_eposID[0] = 2130;
        Original_eposID[1] = 2330;
        Original_pdgID[0] = 3122;
        Original_pdgID[1] = 3312;

        DaughterMassP1[0] = -1;
        DaughterMassP1[1] = -1;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = -1;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }
    else if(DataSetDescr=="LamReso_Xim"){
        pdgID[0] = 0;
        pdgID[1] = 3312;
        Original_eposID[0] = 2130;
        Original_eposID[1] = 2330;
        Original_pdgID[0] = 3122;
        Original_pdgID[1] = 3312;

        DaughterMassP1[0] = Mass_L*0.001;
        DaughterMassP1[1] = Mass_pic*0.001;
        DaughterMassP1[2] = -1;
        DaughterMassP1[3] = -1;

        DaughterMassP2[0] = -1;
        DaughterMassP2[1] = -1;
        DaughterMassP2[2] = -1;
        DaughterMassP2[3] = -1;

        ResoMass[0] = 1.46;
        ResoMass[1] = -1;

        FractionsNbody[0] = 1.0;
        FractionsNbody[1] = 0.0;
        FractionsNbody[2] = 0.0;

        UseResoInfo = true;
    }

    FractionsNbodyC[0] = FractionsNbody[0];
    FractionsNbodyC[1] = FractionsNbodyC[0]+FractionsNbody[1];
    FractionsNbodyC[2] = FractionsNbodyC[1]+FractionsNbody[2];
    if(FractionsNbodyC[2]!=1){
        printf("\033[1;33mWARNING:\033[0m The fractions sum up to %f\n",FractionsNbodyC[2]);
    }

    unsigned NumTotalPairs=0;
    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;
/*
    //angle between (r,k) of the primordials, as a function of k* of the daughters (in CM of the daughters), THE MOST USEFUL I WOULD SAY
    TH2F* h_rkAngle_k_Reso = new TH2F("h_rkAngle_k_Reso","h_rkAngle_k_Reso",64,0,Pi,64,0,1024);
    //angle between (r,k) of the primordials, as a function of k* of the daughters in the rest frame of one of the primordials
    TH2F* h_rkAngle_Mom2 = new TH2F("h_rkAngle_Mom2","h_rkAngle_Mom2",64,0,Pi,64,0,1024);
    //angle between (r,k) of the primordials, as a function of k between them (p_res) of the daughters in the rest frame of one of the primordials
    TH2F* h_rkAngle_MomReso2 = new TH2F("h_rkAngle_MomReso2","h_rkAngle_MomReso2",64,0,Pi,64,0,1024);

    TH2F* h_rkAngle_MomResoA = new TH2F("h_rkAngle_MomResoA","h_rkAngle_MomResoA",64,0,Pi,64,0,1024);
    TH2F* h_rkAngle_MomResoB = new TH2F("h_rkAngle_MomResoB","h_rkAngle_MomResoB",64,0,Pi,64,0,1024);

    //angle between (r*,k*) of the daughters, as a function of k* of the daughters (in CM of the daughters)
    TH2F* h_rkAngle_k_Daughters = new TH2F("h_rkAngle_k_Daughters","h_rkAngle_k_Daughters",64,0,Pi,64,0,1024);
    //angle between (r,k) of the primordials, as a function of k of the primordials (in CM of the primordial)
    TH2F* h_rkAngle_k_Prim = new TH2F("h_rkAngle_k_Prim","h_rkAngle_k_Prim",64,0,Pi,64,0,1024);

    //actually decays the resonances into proton pion or Lambda pion and saves the angles for all initial primordial particles,
    //which after the decay resulted in a pair with k* < FixedValue (say 100 MeV or so to do femto)
    TH1F* h_rkAngle_Reso = new TH1F("h_rkAngle_Reso","h_rkAngle_Reso",64,0,Pi);
*/
    //angle between the parents in the CM of the daughters (function of k*)
    TH2F* h_AP_D = new TH2F("h_AP_D","h_AP_D",64,0,Pi,64,0,1024*1);
    TH2F* h_AP_D_Large_AP_P1 = new TH2F("h_AP_D_Large_AP_P1","h_AP_D_Large_AP_P1",64,0,Pi,64,0,1024*1);
    TH2F* h_AP_D_Large_AP_P2 = new TH2F("h_AP_D_Large_AP_P2","h_AP_D_Large_AP_P2",64,0,Pi,64,0,1024*1);
    TH2F* h_AP_D_LargeAngles = new TH2F("h_AP_D_LargeAngles","h_AP_D_LargeAngles",64,0,Pi,64,0,1024*1);
    TH1F* h1_AP_D = new TH1F("h1_AP_D","h1_AP_D",256,0,Pi);
    //angle between the daughters in the CM of the daughters (function of k*)
    TH2F* h_AD_D = new TH2F("h_AD_D","h_AD_D",64,0,Pi,64,0,1024*1);

    //angle between the parents in the rest frame of the first daughter (function of p_res)
    TH2F* h_AP_D1 = new TH2F("h_AP_D1","h_AP_D1",64,0,Pi,64,0,1024*1);
    //angle between the daughters in the rest frame of the first daughter (function of p_res)
    TH2F* h_AD_D1 = new TH2F("h_AD_D1","h_AD_D1",64,0,Pi,64,0,1024*1);

    //angle between the parents in the rest frame of the second daughter (function of p_res)
    TH2F* h_AP_D2 = new TH2F("h_AP_D2","h_AP_D2",64,0,Pi,64,0,1024*1);
    //angle between the daughters in the rest frame of the second daughter (function of p_res)
    TH2F* h_AD_D2 = new TH2F("h_AD_D2","h_AD_D2",64,0,Pi,64,0,1024*1);

    //angle between the parents in the rest frame of the first parent (function of p_res)
    TH2F* h_AP_P1 = new TH2F("h_AP_P1","h_AP_P1",64,0,Pi,64,0,1024*1);
    TH2F* h_AP_P1_AR_P1P2 = new TH2F("h_AP_P1_AR_P1P2","h_AP_P1_AR_P1P2",64,0,Pi,64,0,Pi);
    TH2F* h_AP_P1_dR_P1P2 = new TH2F("h_AP_P1_dR_P1P2","h_AP_P1_dR_P1P2",64,0,Pi,64,0,64);
    TH2F* h_AP_P1_LargeAngles = new TH2F("h_AP_P1_LargeAngles","h_AP_P1_LargeAngles",64,0,Pi,64,0,1024*1);
    //angle between the daughters in the rest frame of the first parent (function of p_res)
    TH2F* h_AD_P1 = new TH2F("h_AD_P1","h_AD_P1",64,0,Pi,64,0,1024*1);

    //angle between the parents in the rest frame of the second parent (function of p_res)
    TH2F* h_AP_P2 = new TH2F("h_AP_P2","h_AP_P2",64,0,Pi,64,0,1024*1);
    //angle between the daughters in the rest frame of the second parent (function of p_res)
    TH2F* h_AD_P2 = new TH2F("h_AD_P2","h_AD_P2",64,0,Pi,64,0,1024*1);

    TH1F* h_MomDist_P1 = new TH1F("h_MomDist_P1","h_MomDist_P1",64,0,4000);
    TH1F* h_MomDist_P2 = new TH1F("h_MomDist_P2","h_MomDist_P2",64,0,4000);
    TH2F* h_MomDist_P1P2 = new TH2F("h_MomDist_P1P2","h_MomDist_P1P2",64,0,4000,64,0,4000);

    TH1F* h_TauDist_P1 = new TH1F("h_TauDist_P1","h_TauDist_P1",64,0,4);
    TH1F* h_TauDist_P2 = new TH1F("h_TauDist_P2","h_TauDist_P2",64,0,4);
    TH2F* h_TauDist_P1P2 = new TH2F("h_TauDist_P1P2","h_TauDist_P1P2",64,0,4,64,0,4);

    //the angle distribution between the parents in the LAB
    TH1F* h_AP_LAB = new TH1F("h_AP_LAB","h_AP_LAB",256,0,Pi);
    TH2F* h_AP_LAB_AR_LAB = new TH2F("h_AP_LAB_AR_LAB","h_AP_LAB_AR_LAB",64*16,0,Pi,64*16,0,Pi);
    TH2F* h_AP_LAB_AR_LAB_Large_AP_P1 = new TH2F("h_AP_LAB_AR_LAB_Large_AP_P1","h_AP_LAB_AR_LAB_Large_AP_P1",64*16,0,Pi,64*16,0,Pi);
    TH2F* h_AP_LAB_dR_LAB = new TH2F("h_AP_LAB_dR_LAB","h_AP_LAB_dR_LAB",64,0,Pi,64,0,64);

    TH1F* h_AngleRCP1_D = new TH1F("h_AngleRCP1_D","h_AngleRCP1_D",256,0,Pi);
    TH1F* h_AngleRCP2_D = new TH1F("h_AngleRCP2_D","h_AngleRCP2_D",256,0,Pi);
    TH1F* h_AngleP1P2_D = new TH1F("h_AngleP1P2_D","h_AngleP1P2_D",256,0,Pi);

    TH2F* h_AngleRCP1_RCP2_D = new TH2F("h_AngleRCP1_RCP2_D","h_AngleRCP1_RCP2_D",64,0,Pi,64,0,Pi);
    TH2F* h_AngleRCP1_P1P2_D = new TH2F("h_AngleRCP1_P1P2_D","h_AngleRCP1_P1P2_D",64,0,Pi,64,0,Pi);
    TH2F* h_AngleRCP2_P1P2_D   = new TH2F("h_AngleRCP2_P1P2_D","h_AngleRCP2_P1P2_D",64,0,Pi,64,0,Pi);

    TH1F* hResoMass = new TH1F("hResoMass","hResoMass",4096,0,4096);
    TH1I* hNbody = new TH1I("hNbody","hNbody",8,-0.5,7.5);

    TFile* fOut = new TFile(OutFileBaseName+".root","recreate");
    //TFile* fOut = new TFile(OutFileBaseName+"_REPICK.root","recreate");
    TNtuple* InfoTuple = new TNtuple("InfoTuple","InfoTuple","k_D:P1:P2:M1:M2:Tau1:Tau2:AngleRcP1:AngleRcP2:AngleP1P2");
    TNtuple* InfoTuple_ClosePairs = new TNtuple("InfoTuple_ClosePairs","InfoTuple_ClosePairs","k_D:P1:P2:M1:M2:Tau1:Tau2:AngleRcP1:AngleRcP2:AngleP1P2");
    const double ClosePairLimit = 300;

    CatsParticle KittyParticle;
    CatsEvent* KittyEvent = new CatsEvent(pdgID[0],pdgID[1]);
    CatsEvent* KittyFilteredEvent = new CatsEvent(pdgID[0],pdgID[1]);

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;
    TRandom3 rangen(11);

    //GenBod generators, one per particle
    //TGenPhaseSpace* eventRan = new TGenPhaseSpace();
    TGenPhaseSpace* eventRan1=NULL;
    TGenPhaseSpace* eventRan2=NULL;
    //parents
    TLorentzVector TLV_P1;
    TLorentzVector TLV_P2;


    TLorentzVector TLV_Daughter1;
    TLorentzVector TLV_Daughter2;
    TLorentzVector* TLVP_Daughter1;
    TLorentzVector* TLVP_Daughter2;

    //set up the decay channels of each particle
    if(DaughterMassP1[0]>=0) eventRan1 = new TGenPhaseSpace();
    if(DaughterMassP2[0]>=0) eventRan2 = new TGenPhaseSpace();

    //
    while(!feof(InFile)){
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
        TotNumEvents++;
        if(NumPartInEvent>int(HighMultLimit)) {RejectedHighMultEvents++;}
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);

//if(KittyParticle.GetMass()>0.780&&KittyParticle.GetMass()<0.785) {printf("Kicked out an omega\n"); continue;}
            //if(DataSetDescr.Contains("Reso")&&KittyParticle.GetPid()!=Original_pdgID[0]&&KittyParticle.GetPid()!=Original_pdgID[1]&&KittyParticle.GetMass()<DaughterMassP1[0]+DaughterMassP1[1])
            //{
                //printf("mass = %f\n",KittyParticle.GetMass());
                //printf("pid = %i\n",KittyParticle.GetPid());
//!set all masses to the avg
//KittyParticle.SetMass(1.36);
//            }

    //int ParticlePID = KittyParticle.GetPid();
    //to boost statistics we take both the particles and antiparticles and treat them the same
    int ParticlePID = abs(KittyParticle.GetPid());

///Check for neutral pions
//if(ParticlePID==111 || ParticlePID==110){
//if(ParticlePID==131)
//printf("PID=%i; M=%f\n",ParticlePID,KittyParticle.GetMass());
//}
    double ResoWidth;
    int ParentPID;
    //this was used for the baryon baryon femto, where each reso can be mapped to a specific daughter
    //this will not work properly for meson meson, this is why we have to switch it off
    //if(UseResoInfo&&(DataSetDescr=="piReso_Sig"||DataSetDescr=="pi_SigReso"||DataSetDescr=="pi_Sig"||DataSetDescr=="piReso_SigReso")) ResoInfo(ParticlePID,ParentPID,ResoWidth,DataSetDescr);
    //else if(UseResoInfo) ResoInfo(ParticlePID,ParentPID,ResoWidth);
    if(UseResoInfo) ResoInfo(ParticlePID,ParentPID,ResoWidth,DataSetDescr);
    else{
        //to boost statistics we take both the particles and antiparticles and treat them the same
        ParticlePID = abs(ParticlePID);
        KittyParticle.SetPid(ParticlePID);
        //ResoInfo should tell is which is the pid of the daughter associated with a resonances.
        //If it is zero, than it implies the current particle is a primary
        if(ParticlePID==Original_pdgID[0]||ParticlePID==Original_eposID[0]||ParticlePID==Original_pdgID[1]||ParticlePID==Original_eposID[1]){
            ParentPID = 0;
            ResoWidth = 0;
        }
        else if(DataSetDescr=="pi_piReso"||DataSetDescr=="piReso_piReso"){
            //this works only for pions at the moment
            ParentPID = 211;
            ResoWidth = 1;
        }
        //works for pKminus at the moment
        else if(DataSetDescr=="p_KaonReso"||DataSetDescr=="pReso_KaonReso"){
/// Пу дееба, това не работи. За пК, ще се наложи да имаш още работа. Трябва някак си да разбереш дали резонанса идва от протон или от Каон,
/// и на тази база да сложиш ParentPID.
/// Няма да стане :( има дупликати, като н1440, за резонансите на протоните и каоните, ще трябва нещо друго да измислищ
        }
    }
//printf("ParticlePID==%i; ParentPID==%i; ResoWidth==%.2f; %s\n",ParticlePID,ParentPID,ResoWidth,DataSetDescr.Data());
    if(ParticlePID==Original_pdgID[0]||ParticlePID==Original_eposID[0]) {KittyParticle.SetWidth(0);KittyParticle.SetPid(Original_pdgID[0]);}
    else if(ParticlePID==Original_pdgID[1]||ParticlePID==Original_eposID[1]) {KittyParticle.SetWidth(0);KittyParticle.SetPid(Original_pdgID[1]);}
    //if it is a reso
    else if(ParentPID&&ResoWidth){
        KittyParticle.SetWidth(ResoWidth);
        if(ParentPID==Original_pdgID[0]||ParentPID==Original_eposID[0]) {KittyParticle.SetPid(0);}
        else if(ParentPID==Original_pdgID[1]||ParentPID==Original_eposID[1]) {KittyParticle.SetPid(1);}
        //unknown particle
        else {printf(" AAA\n"); continue;}
    }
    //unknown particle
    else {printf(" BBB\n"); continue;}

    //modify the mass of neutrons, which are used as dummies to model f0 and a0 (M = 1 GeV) decays into Kaons
    if(DataSetDescr=="p_KaonReso"||DataSetDescr=="pReso_KaonReso"){
      if(ParticlePID==1220||ParticlePID==2112){
        KittyParticle.SetMass(1.0);
      }
    }

//printf("----\n");
//printf(" id = %i\n",KittyParticle.GetPid());
//printf(" width = %f\n",KittyParticle.GetWidth());

//if(DataSetDescr.Contains("Reso")&&KittyParticle.GetPid()!=Original_pdgID[0]&&KittyParticle.GetPid()!=Original_pdgID[1]
//   &&KittyParticle.GetPid()!=Original_eposID[0]&&KittyParticle.GetPid()!=Original_eposID[1])KittyParticle.SetPid(0);
//else if(KittyParticle.GetPid()==Original_eposID[0])KittyParticle.SetPid(Original_pdgID[0]);
//else if(KittyParticle.GetPid()==Original_eposID[1])KittyParticle.SetPid(Original_pdgID[1]);

//if(DataSetDescr.Contains("Reso")&&KittyParticle.GetPid()!=1120&&
//   KittyParticle.GetPid()!=1114&&KittyParticle.GetPid()!=1128&&KittyParticle.GetPid()!=1228&&KittyParticle.GetPid()!=1127&&KittyParticle.GetPid()!=1227)continue;
//if(DataSetDescr.Contains("Reso")&&KittyParticle.GetPid()!=1120) KittyParticle.SetPid(0);

//printf("pdgID=%i\n",KittyParticle.GetPid());
            //printf("KittyParticle = %f (%.3f,%.3f,%.3f)\n",KittyParticle.GetP(),KittyParticle.GetPx(),KittyParticle.GetPy(),KittyParticle.GetPz());

            //continue;
            //ALICE ACCEPTANCE
            //if(KittyParticle.GetP()<0.4) continue;
            //if(fabs(KittyParticle.GetPseudoRap())>0.8) continue;//!

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;
            if(NumPartInEvent>int(HighMultLimit)) continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            //exclude particle with weird momenta
            if(KittyParticle.GetP()>100){
                continue;
            }


            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type


/*
            //cut out entries at very large rapidity
            if(fabs(KittyParticle.GetPseudoRap())>6){
                continue;
            }
            //cut out pairs that are very far away from the collision point
            if(fabs(KittyParticle.GetR())>16){
                continue;
            }
*/
            KittyEvent->AddParticle(KittyParticle);


if(KittyParticle.GetPid()==0)
hResoMass->Fill(KittyParticle.GetMass()*1000.);
        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)


        for(unsigned uPart=0; uPart<KittyEvent->GetNumParticles1(); uPart++){
            //if(KittyEvent[uMeth][uMix]->GetParticleType1(uPart).GetP()<0.4) continue;
            //if(fabs(KittyEvent[uMeth][uMix]->GetParticleType1(uPart).GetPseudoRap())>0.8) continue;
            KittyFilteredEvent->AddParticle(KittyEvent->GetParticleType1(uPart));
        }
        for(unsigned uPart=0; uPart<KittyEvent->GetNumParticles2(); uPart++){
            //avoid double counting for identical particles
            if(KittyEvent->GetSameType()) break;
//printf("HELLO\n");
            //if(KittyEvent[uMeth][uMix]->GetParticleType2(uPart).GetP()<0.4) continue;
            //if(fabs(KittyEvent[uMeth][uMix]->GetParticleType2(uPart).GetPseudoRap())>0.8) continue;
            KittyFilteredEvent->AddParticle(KittyEvent->GetParticleType2(uPart));
        }
        KittyFilteredEvent->ComputeParticlePairs(false,false);

        NumTotalPairs+=KittyFilteredEvent->GetNumPairs();

        for(unsigned uPair=0; uPair<KittyFilteredEvent->GetNumPairs(); uPair++){
            int pid1 = KittyFilteredEvent->GetParticlePair(uPair).GetParticle(0).GetPid();
            int pid2 = KittyFilteredEvent->GetParticlePair(uPair).GetParticle(1).GetPid();
//printf("pid1/2 = %i/%i\n",pid1,pid2);
            //make sure we have particles of the desired PID
            if((pdgID[0]!=pid1||pdgID[1]!=pid2) && (pdgID[0]!=pid2||pdgID[1]!=pid1)) continue;

            CatsParticle ParticleOriginal1 = KittyFilteredEvent->GetParticlePair(uPair).GetParticle(0);
            CatsParticle ParticleOriginal2 = KittyFilteredEvent->GetParticlePair(uPair).GetParticle(1);

//printf("ParticleOriginal1.GetWidth()=%f(%f)\n",ParticleOriginal1.GetWidth(),KittyFilteredEvent->GetParticlePair(uPair).GetParticle(0).GetWidth());
//printf("ParticleOriginal2.GetWidth()=%f(%f)\n",ParticleOriginal2.GetWidth(),KittyFilteredEvent->GetParticlePair(uPair).GetParticle(1).GetWidth());
//if(ParticleOriginal2.GetMass()<0.5){
//printf("P2 mass = %.2f\n",ParticleOriginal2.GetMass());
//}
            ///correct the mass
            if(pid1<10&&ParticleOriginal1.GetMass()<DaughterMassP1[0]+DaughterMassP1[1]){
//printf("ParticleOriginal1.GetMass() = %f\n",ParticleOriginal1.GetMass());
                ParticleOriginal1.SetMass(ResoMass[0]);
                ParticleOriginal1.Set(ParticleOriginal1.GetT(),ParticleOriginal1.GetX(),ParticleOriginal1.GetY(),ParticleOriginal1.GetZ(),
                                      sqrt(ResoMass[0]*ResoMass[0]+ParticleOriginal1.GetP2()),
                                      ParticleOriginal1.GetPx(),ParticleOriginal1.GetPy(),ParticleOriginal1.GetPz());
            }
            if(pid2<10&&ParticleOriginal2.GetMass()<DaughterMassP2[0]+DaughterMassP2[1]){
//printf("ParticleOriginal2.GetMass() = %f\n",ParticleOriginal2.GetMass());
//printf("ParticleOriginal2.GetMass() = %f\n",ParticleOriginal2.GetMass());
                ParticleOriginal2.SetMass(ResoMass[1]);
                ParticleOriginal2.Set(ParticleOriginal2.GetT(),ParticleOriginal2.GetX(),ParticleOriginal2.GetY(),ParticleOriginal2.GetZ(),
                                      sqrt(ResoMass[1]*ResoMass[1]+ParticleOriginal2.GetP2()),
                                      ParticleOriginal2.GetPx(),ParticleOriginal2.GetPy(),ParticleOriginal2.GetPz());
//printf(" -> new mass: %.2f (%.2f)\n",ResoMass[1],ParticleOriginal2.GetMass());
            }


            CatsParticle Particle1 = ParticleOriginal1;
            CatsParticle Particle2 = ParticleOriginal2;

            h_MomDist_P1->Fill(Particle1.GetP()*1000.);
            h_MomDist_P2->Fill(Particle2.GetP()*1000.);
            h_MomDist_P1P2->Fill(Particle1.GetP()*1000.,Particle2.GetP()*1000.);

            h_TauDist_P1->Fill(Particle1.GetT());
            h_TauDist_P2->Fill(Particle2.GetT());
            h_TauDist_P1P2->Fill(Particle1.GetT(),Particle2.GetT());

            CatsParticle Part1Part2;
            Part1Part2 = Particle1+Particle2;
            CatsLorentzVector RelPart1Part2 = Particle2-Particle1;

        //printf("\nParticle1 = %f (%.3f,%.3f,%.3f) pid %i\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz(),Particle1.GetPid());
        //printf("Particle2 = %f (%.3f,%.3f,%.3f) pid %i\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz(),Particle2.GetPid());

            //the angle (r,k) before any boosts (so in LAB)
            TVector3 vecR;
            TVector3 vecK;
            TVector3 vecP1;
            TVector3 vecP2;
            TVector3 vecRP1;
            TVector3 vecRP2;
            vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());
            double Angle_Lab = vecR.Angle(vecK);
            if(Angle_Lab<0) Angle_Lab += Pi;
            if(Angle_Lab>Pi) Angle_Lab -= Pi;


            //the angle between the two primary particles:
            // in the CM of the daughters
            double Angle_CmDaughters;
            // in the CM of the parents
            double Angle_CmParents;
            // in the rest frame of parent number 1
            double Angle_Parent1;
            // in the rest frame of parent number 2
            double Angle_Parent2;
            // in the rest frame of daughter number 1
            double Angle_Daughter1;
            // in the rest frame of daughter number 2
            double Angle_Daughter2;

            //the same as before, but now the angle is between the two daughters
            double AngleD_CmDaughters;
            // in the CM of the parents
            double AngleD_CmParents;
            // in the rest frame of parent number 1
            double AngleD_Parent1;
            // in the rest frame of parent number 2
            double AngleD_Parent2;
            // in the rest frame of daughter number 1
            double AngleD_Daughter1;
            // in the rest frame of daughter number 2
            double AngleD_Daughter2;

            //the relative momenta between the two primary particles
            // in the CM of the daughters
            double RelMom_CmDaughters;
            // in the CM of the parents
            double RelMom_CmParents;
            // in the rest frame of parent number 1
            double RelMom_Parent1;
            // in the rest frame of parent number 2
            double RelMom_Parent2;
            // in the rest frame of daughter number 1
            double RelMom_Daughter1;
            // in the rest frame of daughter number 2
            double RelMom_Daughter2;

            //the same as before, but now the relative momenta is between the two daughters
            double RelMomD_CmDaughters;
            // in the CM of the parents
            double RelMomD_CmParents;
            // in the rest frame of parent number 1
            double RelMomD_Parent1;
            // in the rest frame of parent number 2
            double RelMomD_Parent2;
            // in the rest frame of daughter number 1
            double RelMomD_Daughter1;
            // in the rest frame of daughter number 2
            double RelMomD_Daughter2;

            //the angle between r_core and the momentum of parent 1 in CM of the daughters
            double AngleRCP1_CmDaughters;
            //the angle between r_core and the momentum of parent 2 in CM of the daughters
            double AngleRCP2_CmDaughters;
            //the angle between the momenta of the parents in CM of the daughters
            double AngleP1P2_CmDaughters;

            double MomP1_CmDaughters;
            double MomP2_CmDaughters;

            TVector3 vP1;
            vP1.SetXYZ(Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            TVector3 vP2;
            vP2.SetXYZ(Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            double Angle_P1P2 = vP1.Angle(vP2);
            //h_kkAngle_Mom2->Fill(Angle_P1P2,(vP1-vP2).Mag()*1000.);
            TVector3 vPosP1;
            vPosP1.SetXYZ(Particle1.GetX(),Particle1.GetY(),Particle1.GetZ());
            TVector3 vPosP2;
            vPosP2.SetXYZ(Particle2.GetX(),Particle2.GetY(),Particle2.GetZ());
            double AngleR_P1P2 = vPosP1.Angle(vPosP2);
            double RelR_P1P2 = (vPosP2-vPosP1).Mag();//r_core

            TLV_P1.SetPxPyPzE(Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz(),Particle1.GetE());
            TLV_P2.SetPxPyPzE(Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz(),Particle2.GetE());

//if(Particle2.GetMass()<0.6||TLV_P2.M()<0.6){
//  double px,py,pz,ener;
//  px = Particle2.GetPx(); py = Particle2.GetPy(); pz = Particle2.GetPz(); ener = Particle2.GetE();
//  printf("M=%.2f (%.2f) [%.2f]\n",Particle2.GetMass(),TLV_P2.M(),sqrt(ener*ener-px*px-py*py-pz*pz));
//}

/*
            if(DaughterMassP1[0]>=0){
                eventRan->SetDecay(TLV_P1, 2, DaughterMassP1);
                eventRan->Generate();
                TLVP_Daughter1 = eventRan->GetDecay(0);
            }
            else{
                TLVP_Daughter1 = &TLV_P1;
            }
            TLV_Daughter1 = *TLVP_Daughter1;

            if(DaughterMassP2[0]>=0){
                eventRan->SetDecay(TLV_P2, 2, DaughterMassP2);
                eventRan->Generate();
                TLVP_Daughter2 = eventRan->GetDecay(0);
            }
            else{
                TLVP_Daughter2 = &TLV_P2;
            }
            TLV_Daughter2 = *TLVP_Daughter2;
*/

            int CrashCounted = 0;

            //decide how many body decay we have
            int NbodyDecayP1=2;
            double NbodyRandomP1;
            double TotalDaughterMassP1 = 1e6;
            //int iterc = 0;
            while(TotalDaughterMassP1>TLV_P1.M()&&TLV_P1.M()>MultiBodyThreshold){
                NbodyRandomP1 = rangen.Uniform();
                for(int inb=2; inb<=4; inb++){
                    if(NbodyRandomP1<=FractionsNbodyC[inb-2]){NbodyDecayP1=inb; break;}
                }
                TotalDaughterMassP1 = DaughterMassP1[0];
                for(int inb=1; inb<NbodyDecayP1; inb++){
                    TotalDaughterMassP1 += DaughterMassP1[inb];
                }
                if(NbodyDecayP1==2 && TotalDaughterMassP1>TLV_P1.M()){
                    printf("WTF!!!!\n");
                    printf(" 1) %f > %f\n",TotalDaughterMassP1,TLV_P1.M());
                }
                //iterc++;
                //if(iterc>100000){
                //  printf("stuck: TLV_P1.M()=%f (%f)\n",TLV_P1.M(),TotalDaughterMassP1);
                //  usleep(250e3);
                //}
            }
            //printf("hi: TLV_P1.M()=%f (%f)\n",TLV_P1.M(),TotalDaughterMassP1);
            //usleep(125e3);

            int NbodyDecayP2=2;
            double NbodyRandomP2;
            double TotalDaughterMassP2 = 1e6;
            while(TotalDaughterMassP2>TLV_P2.M()&&TLV_P2.M()>MultiBodyThreshold){
                NbodyRandomP2 = rangen.Uniform();
                for(int inb=2; inb<=4; inb++){
                    if(NbodyRandomP2<=FractionsNbodyC[inb-2]){NbodyDecayP2=inb; break;}
                }
                TotalDaughterMassP2 = DaughterMassP2[0];
                for(int inb=1; inb<NbodyDecayP2; inb++){
                    TotalDaughterMassP2 += DaughterMassP2[inb];
                }
                if(NbodyDecayP2==2 && TotalDaughterMassP2>TLV_P2.M()){
                    printf("WTF!!!!\n");
                    printf(" 2) %f > %f\n",TotalDaughterMassP2,TLV_P2.M());
                    usleep(1000e3);
                }
            }

            //this could only happen for 2 body decays, when the condition was forced by MultiBodyThreshold
            if(TotalDaughterMassP1==1e6){
              TotalDaughterMassP1 = DaughterMassP1[0]+DaughterMassP1[1];
            }
            if(TotalDaughterMassP2==1e6){
              TotalDaughterMassP2 = DaughterMassP2[0]+DaughterMassP2[1];
            }

            hNbody->Fill(NbodyDecayP1);
            hNbody->Fill(NbodyDecayP2);

            if(eventRan1) eventRan1->SetDecay(TLV_P1, NbodyDecayP1, DaughterMassP1);
            if(eventRan2) eventRan2->SetDecay(TLV_P2, NbodyDecayP2, DaughterMassP2);

            //if(eventRan1) printf("1\n");
            //if(eventRan2) printf("2\n");

            const unsigned NumRepicks = PickBestDecay?16:1;
            TLorentzVector* Sol_D1 = new TLorentzVector [NumRepicks];
            TLorentzVector* Sol_D2 = new TLorentzVector [NumRepicks];
            unsigned RepickCounter = 0;
            unsigned BestSol = 0;
            double BestRelMom = 1e6;

            //if(PickBestDecay){
                for(unsigned uRepick=0; uRepick<NumRepicks; uRepick++){
                    //the resonances are decayed here
                    if(eventRan1) eventRan1->Generate();
                    if(eventRan2) eventRan2->Generate();

                    //get the daughters
                    if(eventRan1) TLVP_Daughter1 = eventRan1->GetDecay(0);
                    else TLVP_Daughter1 = &TLV_P1;
                    if(eventRan2) TLVP_Daughter2 = eventRan2->GetDecay(0);
                    else TLVP_Daughter2 = &TLV_P2;

                    Sol_D1[RepickCounter] = *TLVP_Daughter1;
                    Sol_D2[RepickCounter] = *TLVP_Daughter2;

                    TLorentzVector D1D2_diff = Sol_D1[RepickCounter]-Sol_D2[RepickCounter];

                    if(D1D2_diff.P()<BestRelMom){
                        BestRelMom = D1D2_diff.P();
                        BestSol = RepickCounter;
                    }

                    RepickCounter++;
                }
            //}

            TLV_Daughter1 = Sol_D1[BestSol];
//printf("TLV_Daughter2 set to %u\n",BestSol);
            TLV_Daughter2 = Sol_D2[BestSol];

            delete [] Sol_D1;
            delete [] Sol_D2;

///////////////////////////////////////////////////////////////////////////
            ////the resonances are decayed here
            //if(eventRan1) eventRan1->Generate();
            //if(eventRan2) eventRan2->Generate();

            if(TLV_P1.M()<=TotalDaughterMassP1){
                //printf("TLV_P1.M()=%f\n",TLV_P1.M());
                //abort();
printf("TROUBLE EXPECTED...\n");
printf("TLV_P1.M()=%f <= %f\n",TLV_P1.M(),TotalDaughterMassP1);
            }
            if(TLV_P2.M()<=TotalDaughterMassP2){
                //printf("TLV_P2.M()=%f\n",TLV_P2.M());
                //abort();
printf("TROUBLE EXPECTED...\n");
printf("TLV_P2.M()=%f <= %f\n",TLV_P2.M(),TotalDaughterMassP2);
            }

            ////get the daughters
            //if(eventRan1) TLVP_Daughter1 = eventRan1->GetDecay(0);
            //else TLVP_Daughter1 = &TLV_P1;
            //if(eventRan2) TLVP_Daughter2 = eventRan2->GetDecay(0);
            //else TLVP_Daughter2 = &TLV_P2;

            //TLV_Daughter1 = *TLVP_Daughter1;
            //TLV_Daughter2 = *TLVP_Daughter2;
//////////////////////////////////////////////////////////////////////////


            //TLV_Daughter1.Print();
            //TLV_Daughter2.Print();
if(TLV_Daughter2.Px()!=TLV_Daughter2.Px()){
printf("TROUBLE INCOMING!\n");
printf("eventRan1=%p\n",eventRan1);
printf("eventRan2=%p\n",eventRan2);
printf("TLV_P2.M()=%f\n",TLV_P2.M());
printf("TotalDaughterMassP2=%f\n",TotalDaughterMassP2);
printf("TLV_Daughter2 set to %u\n",BestSol);
printf("NbodyDecayP2 set to %u\n",NbodyDecayP2);
}
            CatsParticle DaughterOriginal1;
            CatsParticle DaughterOriginal2;
            DaughterOriginal1.Set(Particle1.GetT(),Particle1.GetX(),Particle1.GetY(),Particle1.GetZ(),TLV_Daughter1.E(),TLV_Daughter1.Px(),TLV_Daughter1.Py(),TLV_Daughter1.Pz());
            DaughterOriginal2.Set(Particle2.GetT(),Particle2.GetX(),Particle2.GetY(),Particle2.GetZ(),TLV_Daughter2.E(),TLV_Daughter2.Px(),TLV_Daughter2.Py(),TLV_Daughter2.Pz());
            CatsParticle Daughter1 = DaughterOriginal1;
            CatsParticle Daughter2 = DaughterOriginal2;

            TLorentzVector TVL_Daughters = TLV_Daughter1 + TLV_Daughter2;

            //!THE NEW STUFF
            //boost into the rest frame of the daughters
            CatsLorentzVector BoostVector;
            BoostVector.Set(0,0,0,0,TVL_Daughters.E(),TVL_Daughters.Px(),TVL_Daughters.Py(),TVL_Daughters.Pz());
            //printf("BV = (%.3f,%.3f,%.3f)\n",BoostVector.GetPx(),BoostVector.GetPy(),BoostVector.GetPz());
            Particle1.Boost(BoostVector);
            Particle2.Boost(BoostVector);
            Daughter1.Boost(BoostVector);
            Daughter2.Boost(BoostVector);
        //printf("\np1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf("p2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
        //printf("d1 = %f (%.3f,%.3f,%.3f)\n",Daughter1.GetP(),Daughter1.GetPx(),Daughter1.GetPy(),Daughter1.GetPz());
        //printf("d2 = %f (%.3f,%.3f,%.3f)\n",Daughter2.GetP(),Daughter2.GetPx(),Daughter2.GetPy(),Daughter2.GetPz());
            vecR.SetXYZ(Particle2.GetX()-Particle1.GetX(),Particle2.GetY()-Particle1.GetY(),Particle2.GetZ()-Particle1.GetZ());//r_core
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());//k* of the parents
            TVector3 vecR_Daughters;
            vecR_Daughters.SetXYZ(Daughter1.GetX()-Daughter2.GetX(),Daughter1.GetY()-Daughter2.GetY(),Daughter1.GetZ()-Daughter2.GetZ());
            TVector3 vecK_Daughters;
            vecK_Daughters.SetXYZ(Daughter2.GetPx()-Daughter1.GetPx(),Daughter2.GetPy()-Daughter1.GetPy(),Daughter2.GetPz()-Daughter1.GetPz());
            vecRP1.SetXYZ(Particle1.GetX(),Particle1.GetY(),Particle1.GetZ());
            vecRP2.SetXYZ(Particle2.GetX(),Particle2.GetY(),Particle2.GetZ());
            vecP1.SetXYZ(Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            vecP2.SetXYZ(Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            Angle_CmDaughters = vecR.Angle(vecK);
            if(Angle_CmDaughters<0) Angle_CmDaughters += Pi;
            if(Angle_CmDaughters>Pi) Angle_CmDaughters -= Pi;

            AngleD_CmDaughters = vecR_Daughters.Angle(vecK_Daughters);
            if(AngleD_CmDaughters<0) AngleD_CmDaughters += Pi;
            if(AngleD_CmDaughters>Pi) AngleD_CmDaughters -= Pi;

            RelMom_CmDaughters = vecK.Mag()*1000.;
            RelMomD_CmDaughters = vecK_Daughters.Mag()*1000.;

            //if the core is defined as r_c = rp2 - rp1, than the the sign convention is as noted on the tablet
            AngleRCP1_CmDaughters = vecR.Angle(vecP1);
            //if(vecR.Dot(vecP1)<0)AngleRCP1_CmDaughters=-AngleRCP1_CmDaughters;
            AngleRCP2_CmDaughters = vecR.Angle(vecP2);
            //if(vecR.Dot(vecP2)<0)AngleRCP2_CmDaughters=-AngleRCP2_CmDaughters;
            AngleP1P2_CmDaughters = vecP1.Angle(vecP2);
            //if(vecP1.Dot(vecP2)<0)AngleP1P2_CmDaughters=-AngleP1P2_CmDaughters;

            TVector3 R_CORE = vecRP2-vecRP1;
            TVector3 BGT1; BGT1.SetXYZ(Particle1.GetPx()/1.400*1.5,Particle1.GetPy()/1.400*1.5,Particle1.GetPz()/1.400*1.5);
            TVector3 BGT2; BGT2.SetXYZ(Particle2.GetPx()/1.400*1.5,Particle2.GetPy()/1.400*1.5,Particle2.GetPz()/1.400*1.5);
            TVector3 R_EFF = R_CORE - BGT1 + BGT2;
            //if(RelMomD_CmDaughters<200){
            //printf("\n");
            //printf("K* = %f\n",RelMomD_CmDaughters);
            //printf("R_CORE = %f\n",R_CORE.Mag());
            //printf("R_EFF = %f\n",R_EFF.Mag());
            //}


//double D_COS = (vecR.X()*vecP1.X()+vecR.Y()*vecP1.Y()+vecR.Z()*vecP1.Z())/(vecR.Mag()*vecP1.Mag());
//double R_COS = cos(AngleRCP1_CmDaughters);
//double D_A = acos(D_COS);
//double R_A = AngleRCP1_CmDaughters;
//if(D_A*RadToDeg>180){
//printf("\n");
//printf("Dimi: A=%.3f; C=%.3f\n",D_A*RadToDeg,D_COS);
//printf("Root: A=%.3f; C=%.3f\n",R_A*RadToDeg,R_COS);
//}




            MomP1_CmDaughters = vecP1.Mag();
            MomP2_CmDaughters = vecP2.Mag();

        //printf("p1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf("p2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
        //printf("Angle = %f\n",Angle*RadToDeg);
        //printf("RelMomD_CmDaughters = %f\n",RelMomD_CmDaughters);

            //angle between (r,k) of the primordials, as a function of k* of the daughters
            //h_rkAngle_k_Reso->Fill(Angle,RelMomD_CmDaughters);
            //h_rkAngle_k_Daughters->Fill(AngleDaughters,RelMomD_CmDaughters);
            //if(RelMomD_CmDaughters<200) h_rkAngle_Reso->Fill(Angle);
            //if(Angle>13.14){
            //    printf("----------------\n");
            //    printf(" Angle = %f\n\n",Angle);
            //    printf(" RelMomD_CmDaughters = %f\n\n",RelMomD_CmDaughters);
            //    printf(" p1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            //    printf(" p2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            //    printf(" o1 = %f (%.3f,%.3f,%.3f)\n",ParticleOriginal1.GetP(),ParticleOriginal1.GetPx(),ParticleOriginal1.GetPy(),ParticleOriginal1.GetPz());
            //    printf(" o2 = %f (%.3f,%.3f,%.3f)\n",ParticleOriginal2.GetP(),ParticleOriginal2.GetPx(),ParticleOriginal2.GetPy(),ParticleOriginal2.GetPz());
            //    printf(" d1 = %f (%.3f,%.3f,%.3f)\n",Daughter1.GetP(),Daughter1.GetPx(),Daughter1.GetPy(),Daughter1.GetPz());
            //    printf(" d2 = %f (%.3f,%.3f,%.3f)\n",Daughter2.GetP(),Daughter2.GetPx(),Daughter2.GetPy(),Daughter2.GetPz());
            //}


            //! get in the CM of the primordials for a change
            Particle1 = ParticleOriginal1;
            Particle2 = ParticleOriginal2;
            Daughter1 = DaughterOriginal1;
            Daughter2 = DaughterOriginal2;
        //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            BoostVector.Set(0,0,0,0,Part1Part2.GetE(),Part1Part2.GetPx(),Part1Part2.GetPy(),Part1Part2.GetPz());
            Particle1.Boost(BoostVector);
            Particle2.Boost(BoostVector);
            Daughter1.Boost(BoostVector);
            Daughter2.Boost(BoostVector);

            vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());
            vecR_Daughters.SetXYZ(Daughter1.GetX()-Daughter2.GetX(),Daughter1.GetY()-Daughter2.GetY(),Daughter1.GetZ()-Daughter2.GetZ());
            vecK_Daughters.SetXYZ(Daughter2.GetPx()-Daughter1.GetPx(),Daughter2.GetPy()-Daughter1.GetPy(),Daughter2.GetPz()-Daughter1.GetPz());
            Angle_CmParents = vecR.Angle(vecK);
            if(Angle_CmParents<0) Angle_CmParents += Pi;
            if(Angle_CmParents>Pi) Angle_CmParents -= Pi;

            AngleD_CmParents = vecR_Daughters.Angle(vecK_Daughters);
            if(AngleD_CmParents<0) AngleD_CmParents += Pi;
            if(AngleD_CmParents>Pi) AngleD_CmParents -= Pi;

            RelMom_CmParents = vecK.Mag()*1000.;
            RelMomD_CmParents = vecK_Daughters.Mag()*1000.;

            //h_rkAngle_k_Prim->Fill(Angle_CmParents,RelMom_CmParents);


            //!Go to the rest frame of ONE of the daughters (1st) and plot the angle between the other guy (p_res) and r
            Particle1 = ParticleOriginal1;
            Particle2 = ParticleOriginal2;
            Daughter1 = DaughterOriginal1;
            Daughter2 = DaughterOriginal2;
            //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            BoostVector.Set(0,0,0,0,Daughter1.GetE(),Daughter1.GetPx(),Daughter1.GetPy(),Daughter1.GetPz());
            Particle1.Boost(BoostVector);
            Particle2.Boost(BoostVector);
            Daughter1.Boost(BoostVector);
            Daughter2.Boost(BoostVector);
        //printf("1st daughter:\n");
        //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());
            vecR_Daughters.SetXYZ(Daughter1.GetX()-Daughter2.GetX(),Daughter1.GetY()-Daughter2.GetY(),Daughter1.GetZ()-Daughter2.GetZ());
            vecK_Daughters.SetXYZ(Daughter2.GetPx()-Daughter1.GetPx(),Daughter2.GetPy()-Daughter1.GetPy(),Daughter2.GetPz()-Daughter1.GetPz());
            Angle_Daughter1 = vecR.Angle(vecK);
            if(Angle_Daughter1<0) Angle_Daughter1 += Pi;
            if(Angle_Daughter1>Pi) Angle_Daughter1 -= Pi;
            RelMom_Daughter1 = vecK.Mag()*1000.;
            RelMomD_Daughter1 = vecK_Daughters.Mag()*1000.;
            //h_rkAngle_MomResoA->Fill(Angle_Daughter1,RelMom_Daughter1);
            //printf("A: A=%.3f; k=%.0f\n",Angle_Daughter1,RelMom_Daughter1);

            //!Go to the rest frame of ONE of the daughters (2nd) and plot the angle between the other guy (p_res) and r
            Particle1 = ParticleOriginal1;
            Particle2 = ParticleOriginal2;
            Daughter1 = DaughterOriginal1;
            Daughter2 = DaughterOriginal2;
            BoostVector.Set(0,0,0,0,Daughter2.GetE(),Daughter2.GetPx(),Daughter2.GetPy(),Daughter2.GetPz());
            Particle1.Boost(BoostVector);
            Particle2.Boost(BoostVector);
            Daughter1.Boost(BoostVector);
            Daughter2.Boost(BoostVector);
        //printf("2nd daughter:\n");
        //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());
            vecR_Daughters.SetXYZ(Daughter1.GetX()-Daughter2.GetX(),Daughter1.GetY()-Daughter2.GetY(),Daughter1.GetZ()-Daughter2.GetZ());
            vecK_Daughters.SetXYZ(Daughter2.GetPx()-Daughter1.GetPx(),Daughter2.GetPy()-Daughter1.GetPy(),Daughter2.GetPz()-Daughter1.GetPz());
            Angle_Daughter2 = vecR.Angle(vecK);
            if(Angle_Daughter2<0) Angle_Daughter2 += Pi;
            if(Angle_Daughter2>Pi) Angle_Daughter2 -= Pi;
            RelMom_Daughter2 = vecK.Mag()*1000.;
            RelMomD_Daughter2 = vecK_Daughters.Mag()*1000.;
            //h_rkAngle_MomResoB->Fill(Angle_Daughter2,RelMom_Daughter2);
            //printf("B: A=%.3f; k=%.0f\n",Angle_Daughter2,RelMom_Daughter2);

            //!Go to the rest frame of ONE of the parents (1st) and plot the angle between the other guy and r
            Particle1 = ParticleOriginal1;
            Particle2 = ParticleOriginal2;
            Daughter1 = DaughterOriginal1;
            Daughter2 = DaughterOriginal2;
            //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            BoostVector.Set(0,0,0,0,Particle1.GetE(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            Particle1.Boost(BoostVector);
            Particle2.Boost(BoostVector);
            Daughter1.Boost(BoostVector);
            Daughter2.Boost(BoostVector);
        //printf("1st daughter:\n");
        //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());
            vecR_Daughters.SetXYZ(Daughter1.GetX()-Daughter2.GetX(),Daughter1.GetY()-Daughter2.GetY(),Daughter1.GetZ()-Daughter2.GetZ());
            vecK_Daughters.SetXYZ(Daughter2.GetPx()-Daughter1.GetPx(),Daughter2.GetPy()-Daughter1.GetPy(),Daughter2.GetPz()-Daughter1.GetPz());
            Angle_Parent1 = vecR.Angle(vecK);
            if(Angle_Parent1<0) Angle_Parent1 += Pi;
            if(Angle_Parent1>Pi) Angle_Parent1 -= Pi;
            RelMom_Parent1 = vecK.Mag()*1000.;
            RelMomD_Parent1 = vecK_Daughters.Mag()*1000.;
            //h_rkAngle_MomResoA->Fill(Angle_Parent1,RelMom_Parent1);
            //printf("A: A=%.3f; k=%.0f\n",Angle_Parent1,RelMom_Parent1);

            //!Go to the rest frame of ONE of the parents (2nd) and plot the angle between the other guy and r
            Particle1 = ParticleOriginal1;
            Particle2 = ParticleOriginal2;
            Daughter1 = DaughterOriginal1;
            Daughter2 = DaughterOriginal2;
            BoostVector.Set(0,0,0,0,Particle2.GetE(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            Particle1.Boost(BoostVector);
            Particle2.Boost(BoostVector);
            Daughter1.Boost(BoostVector);
            Daughter2.Boost(BoostVector);
        //printf("2nd daughter:\n");
        //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());
            vecR_Daughters.SetXYZ(Daughter1.GetX()-Daughter2.GetX(),Daughter1.GetY()-Daughter2.GetY(),Daughter1.GetZ()-Daughter2.GetZ());
            vecK_Daughters.SetXYZ(Daughter2.GetPx()-Daughter1.GetPx(),Daughter2.GetPy()-Daughter1.GetPy(),Daughter2.GetPz()-Daughter1.GetPz());
            Angle_Parent2 = vecR.Angle(vecK);
            if(Angle_Parent2<0) Angle_Parent2 += Pi;
            if(Angle_Parent2>Pi) Angle_Parent2 -= Pi;
            RelMom_Parent2 = vecK.Mag()*1000.;
            RelMomD_Parent2 = vecK_Daughters.Mag()*1000.;
            //h_rkAngle_MomResoB->Fill(Angle_Parent2,RelMom_Parent2);
            //printf("B: A=%.3f; k=%.0f\n",Angle_Parent2,RelMom_Parent2);

/*
            //!the old stuff for PRIM-RESO
            //BoostVector.Set(0,0,0,0,BoostVector.GetE(),-BoostVector.GetPx(),-BoostVector.GetPy(),-BoostVector.GetPz());
            Particle1 = ParticleOriginal1;
            Particle2 = ParticleOriginal2;
            Daughter1 = DaughterOriginal1;
            Daughter2 = DaughterOriginal2;
        //printf("\np1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf("p2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
        //printf("d1 = %f (%.3f,%.3f,%.3f)\n",Daughter1.GetP(),Daughter1.GetPx(),Daughter1.GetPy(),Daughter1.GetPz());
        //printf("d2 = %f (%.3f,%.3f,%.3f)\n",Daughter2.GetP(),Daughter2.GetPx(),Daughter2.GetPy(),Daughter2.GetPz());
            if(DataSetDescr=="p_p"){
                BoostVector.Set(0,0,0,0,0,0,0,0);
            }
            else if(DataSetDescr=="p_pReso"){
                BoostVector.Set(0,0,0,0,Particle1.GetE(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            }
            else if(DataSetDescr=="pReso_pReso"){
                BoostVector.Set(0,0,0,0,Daughter1.GetE(),Daughter1.GetPx(),Daughter1.GetPy(),Daughter1.GetPz());
                //BoostVector.Set(0,0,0,0,TVL_Daughters.E(),TVL_Daughters.Px(),TVL_Daughters.Py(),TVL_Daughters.Pz());
            }
            //p-Lambda correlation
            else if(DataSetDescr=="p_LamReso"){
                BoostVector.Set(0,0,0,0,Particle1.GetE(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
            }
            else if(DataSetDescr=="pReso_Lam"){
                BoostVector.Set(0,0,0,0,Particle2.GetE(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            }
            else if(DataSetDescr=="pReso_LamReso"){
                BoostVector.Set(0,0,0,0,Daughter1.GetE(),Daughter1.GetPx(),Daughter1.GetPy(),Daughter1.GetPz());
            }
            Particle1.Boost(BoostVector);
            Particle2.Boost(BoostVector);
            Daughter1.Boost(BoostVector);
            Daughter2.Boost(BoostVector);
        //printf("PRIM-RESO:\n");
        //printf(" Particle1 = %f (%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
        //printf(" Particle2 = %f (%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
            vecR.SetXYZ(Particle1.GetX()-Particle2.GetX(),Particle1.GetY()-Particle2.GetY(),Particle1.GetZ()-Particle2.GetZ());
            vecK.SetXYZ(Particle2.GetPx()-Particle1.GetPx(),Particle2.GetPy()-Particle1.GetPy(),Particle2.GetPz()-Particle1.GetPz());
            vecK_Daughters.SetXYZ(Daughter2.GetPx()-Daughter1.GetPx(),Daughter2.GetPy()-Daughter1.GetPy(),Daughter2.GetPz()-Daughter1.GetPz());
            Angle = vecR.Angle(vecK);
            if(Angle<0) Angle += Pi;
            if(Angle>Pi) Angle -= Pi;
            relK = vecK.Mag()*1000.;
            relK_Daughters = vecK_Daughters.Mag()*1000.;
            //angle between (r,k) of the primordials, as a function of k* of the daughters in the rest frame of one of the primordials
            h_rkAngle_Mom2->Fill(Angle,relK_Daughters);
            //angle between (r,k) of the primordials, as a function of k between them (p_res) of the daughters in the rest frame of one of the primordials
            //if(Angle<3.1)
            h_rkAngle_MomReso2->Fill(Angle,relK);
//if(Angle<0.001){
if(Angle>13.1415){
printf("\np1 = %f(%.3f,%.3f,%.3f) r %f(%.3f,%.3f,%.3f)\n",Particle1.GetP(),Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz(),
       Particle1.GetR(),Particle1.GetX(),Particle1.GetY(),Particle1.GetZ());
printf("p2 = %f(%.3f,%.3f,%.3f) r %f(%.3f,%.3f,%.3f)\n",Particle2.GetP(),Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz(),
       Particle2.GetR(),Particle2.GetX(),Particle2.GetY(),Particle2.GetZ());
printf("d1 = %f(%.3f,%.3f,%.3f) r %f(%.3f,%.3f,%.3f)\n",Daughter1.GetP(),Daughter1.GetPx(),Daughter1.GetPy(),Daughter1.GetPz(),
       Daughter1.GetR(),Daughter1.GetX(),Daughter1.GetY(),Daughter1.GetZ());
printf("d2 = %f(%.3f,%.3f,%.3f) r %f(%.3f,%.3f,%.3f)\n",Daughter2.GetP(),Daughter2.GetPx(),Daughter2.GetPy(),Daughter2.GetPz(),
       Daughter2.GetR(),Daughter2.GetX(),Daughter2.GetY(),Daughter2.GetZ());


TVector3 UnitPart1;
double MAG1 = ParticleOriginal1.GetP();
UnitPart1.SetXYZ(ParticleOriginal1.GetPx()/MAG1,ParticleOriginal1.GetPy()/MAG1,ParticleOriginal1.GetPz()/MAG1);
TVector3 UnitPart2;
double MAG2 = ParticleOriginal2.GetP();
UnitPart2.SetXYZ(ParticleOriginal2.GetPx()/MAG2,ParticleOriginal2.GetPy()/MAG2,ParticleOriginal2.GetPz()/MAG2);

printf("IN LAB:\np1 = (%.3f,%.3f,%.3f)\n",UnitPart1.X(),UnitPart1.Y(),UnitPart1.Z());
printf("p2 = (%.3f,%.3f,%.3f)\n",UnitPart2.X(),UnitPart2.Y(),UnitPart2.Z());
TVector3 vecOriginalR;
vecOriginalR.SetXYZ(ParticleOriginal1.GetX()-ParticleOriginal2.GetX(),ParticleOriginal1.GetY()-ParticleOriginal2.GetY(),ParticleOriginal1.GetZ()-ParticleOriginal2.GetZ());
printf("deltaR = %.3f\n",vecOriginalR.Mag());
}
*/

            if(Angle_Parent1>13.14159){
                printf("Back to back pairs:\n");
printf("p1 = %f(%.3f,%.3f) r %f(%.3f,%.3f)\n",ParticleOriginal1.GetP(),ParticleOriginal1.GetPtheta(),ParticleOriginal1.GetPphi(),
       ParticleOriginal1.GetR(),ParticleOriginal1.GetTheta(),ParticleOriginal1.GetPhi());
printf("p2 = %f(%.3f,%.3f) r %f(%.3f,%.3f)\n",ParticleOriginal2.GetP(),ParticleOriginal2.GetPtheta(),ParticleOriginal2.GetPphi(),
       ParticleOriginal2.GetR(),ParticleOriginal2.GetTheta(),ParticleOriginal2.GetPhi());
                printf(" delta_P_Parents = %.3f\n",RelPart1Part2.GetP());
                printf(" delta_Pt_Parents = %.3f\n",RelPart1Part2.GetPt());

            }


           // if(Angle_P1P2<0.02&&AngleR_P1P2<0.02)
            h_AP_D->Fill(Angle_CmDaughters,RelMomD_CmDaughters);
            if(Angle_P1P2>1.57&&AngleR_P1P2>1.57)
            h_AP_D_LargeAngles->Fill(Angle_CmDaughters,RelMomD_CmDaughters);
            if(RelMomD_CmDaughters<200) h1_AP_D->Fill(Angle_CmDaughters);
            if(Angle_Parent1>3.138) h_AP_D_Large_AP_P1->Fill(Angle_CmDaughters,RelMomD_CmDaughters);
            if(Angle_Parent2>3.138) h_AP_D_Large_AP_P2->Fill(Angle_CmDaughters,RelMomD_CmDaughters);
            h_AD_D->Fill(AngleD_CmDaughters,RelMomD_CmDaughters);

            h_AP_D1->Fill(Angle_Daughter1,RelMomD_Daughter1);
            h_AD_D1->Fill(AngleD_Daughter1,RelMomD_Daughter1);

            h_AP_D2->Fill(Angle_Daughter2,RelMomD_Daughter2);
            h_AD_D2->Fill(AngleD_Daughter2,RelMomD_Daughter2);

            //if(RelPart1Part2.GetP()>20&&RelPart1Part2.GetPt()<0.5)
            //if(Angle_P1P2<0.02&&AngleR_P1P2<0.02)
            h_AP_P1->Fill(Angle_Parent1,RelMomD_Parent1);
            if(Angle_P1P2>1.57&&AngleR_P1P2>1.57)
            h_AP_P1_LargeAngles->Fill(Angle_Parent1,RelMomD_Parent1);
            h_AP_P1_AR_P1P2->Fill(Angle_Parent1,AngleR_P1P2);
            h_AP_P1_dR_P1P2->Fill(Angle_Parent1,RelR_P1P2);

            h_AD_P1->Fill(AngleD_Parent1,RelMomD_Parent1);

            h_AP_P2->Fill(Angle_Parent2,RelMomD_Parent2);
            h_AD_P2->Fill(AngleD_Parent2,RelMomD_Parent2);



            //if(Angle_Parent1>3.1&&RelMomD_Parent1<500)
            h_AP_LAB->Fill(Angle_P1P2);
            h_AP_LAB_AR_LAB->Fill(Angle_P1P2,AngleR_P1P2);
            if(Angle_Parent1>3.1&&RelMomD_Parent1<1000)
            h_AP_LAB_AR_LAB_Large_AP_P1->Fill(Angle_P1P2,AngleR_P1P2);
            //if(Angle_Parent1>3.1&&RelMomD_Parent1<1000)
            h_AP_LAB_dR_LAB->Fill(Angle_P1P2,RelR_P1P2);

            h_AngleRCP1_D->Fill(AngleRCP1_CmDaughters);
            h_AngleRCP2_D->Fill(AngleRCP2_CmDaughters);
            h_AngleP1P2_D->Fill(AngleP1P2_CmDaughters);
            h_AngleRCP1_RCP2_D->Fill(AngleRCP1_CmDaughters,AngleRCP2_CmDaughters);
            h_AngleRCP1_P1P2_D->Fill(AngleRCP1_CmDaughters,AngleP1P2_CmDaughters);
            h_AngleRCP2_P1P2_D->Fill(AngleRCP2_CmDaughters,AngleP1P2_CmDaughters);

            if(Angle_Parent1>3.1415&&Angle_P1P2>3.0&&AngleR_P1P2>3.0){
                //printf("SUPER-BACK-TO-BACK:\n");
                //printf( "Parents k* = %.1f\n",RelMom_CmParents);
                //printf( "Daughters k* = %.1f\n",RelMomD_CmDaughters);
            }
            Float_t CONTAINER[10];
            CONTAINER[0] = RelMomD_CmDaughters;
            CONTAINER[1] = MomP1_CmDaughters*1000.;
            CONTAINER[2] = MomP2_CmDaughters*1000.;
            CONTAINER[3] = ParticleOriginal1.GetMass()*1000.;
            CONTAINER[4] = ParticleOriginal2.GetMass()*1000.;
            CONTAINER[5] = ParticleOriginal1.GetWidth()?NuToFm/(ParticleOriginal1.GetWidth()*1000.):0;//in fm
            CONTAINER[6] = ParticleOriginal2.GetWidth()?NuToFm/(ParticleOriginal2.GetWidth()*1000.):0;//in fm
//printf("ParticleOriginal1.GetWidth()=%f\n",ParticleOriginal1.GetWidth());
//printf("ParticleOriginal2.GetWidth()=%f\n",ParticleOriginal2.GetWidth());
            CONTAINER[7] = AngleRCP1_CmDaughters;
            CONTAINER[8] = AngleRCP2_CmDaughters;
            CONTAINER[9] = AngleP1P2_CmDaughters;


const double r_core = R_CORE.Mag();
const double bgt1 = MomP1_CmDaughters/1.4*1.5;
const double bgt2 = MomP2_CmDaughters/1.4*1.5;
const double cosR1 = cos(AngleRCP1_CmDaughters);
const double cosR2 = cos(AngleRCP2_CmDaughters);
const double cos12 = cos(AngleP1P2_CmDaughters);
//if(RelMomD_CmDaughters<200){
//printf("The output from the container:\n");
//printf(" k* = %f\n",RelMomD_CmDaughters);
//printf(" r_core = %f\n",r_core);
//printf(" r^2_eff = %f\n",r_core*r_core+bgt1*bgt1+bgt2*bgt2-2.*r_core*bgt1*cosR1+2.*r_core*bgt2*cosR2-2.*bgt1*bgt2*cos12);
//printf(" r_eff = %f\n",sqrt(r_core*r_core+bgt1*bgt1+bgt2*bgt2-2.*r_core*bgt1*cosR1+2.*r_core*bgt2*cosR2-2.*bgt1*bgt2*cos12));
//}
if(r_core*r_core+bgt1*bgt1+bgt2*bgt2-2.*r_core*bgt1*cosR1+2.*r_core*bgt2*cosR2-2.*bgt1*bgt2*cos12<0){
    printf("PROBLEM\n");
}
//! THE FIRST PARTICLE IS A PROTON IN CASE OF pReso
//printf("%.3f %.3f %.3f %.3f\n",CONTAINER[0],CONTAINER[1],CONTAINER[2],CONTAINER[3]);
//printf("pid1=%i\n",pid1);
//printf("pid2=%i\n",pid2);


if(RelMomD_CmDaughters!=RelMomD_CmDaughters || AngleRCP1_CmDaughters!=AngleRCP1_CmDaughters ||
   AngleRCP2_CmDaughters!=AngleRCP2_CmDaughters || AngleP1P2_CmDaughters!=AngleP1P2_CmDaughters){
printf("\nRelMom_CmDaughters=%f\n",RelMomD_CmDaughters);
printf("AngleRCP1_CmDaughters=%f\n",AngleRCP1_CmDaughters);
printf("AngleRCP2_CmDaughters=%f\n",AngleRCP2_CmDaughters);
printf("AngleP1P2_CmDaughters=%f\n",AngleP1P2_CmDaughters);
printf("p1 = %f(%.3f,%.3f) r %f(%.3f,%.3f)\n",ParticleOriginal1.GetP(),ParticleOriginal1.GetPtheta(),ParticleOriginal1.GetPphi(),
       ParticleOriginal1.GetR(),ParticleOriginal1.GetTheta(),ParticleOriginal1.GetPhi());
printf("p2 = %f(%.3f,%.3f) r %f(%.3f,%.3f)\n",ParticleOriginal2.GetP(),ParticleOriginal2.GetPtheta(),ParticleOriginal2.GetPphi(),
       ParticleOriginal2.GetR(),ParticleOriginal2.GetTheta(),ParticleOriginal2.GetPhi());
printf("d1 = %f(%.3f,%.3f) r %f(%.3f,%.3f)\n",DaughterOriginal1.GetP(),DaughterOriginal1.GetPtheta(),DaughterOriginal1.GetPphi(),
       DaughterOriginal1.GetR(),DaughterOriginal1.GetTheta(),DaughterOriginal1.GetPhi());
printf("d2 = %f(%.3f,%.3f) r %f(%.3f,%.3f)\n",DaughterOriginal2.GetP(),DaughterOriginal2.GetPtheta(),DaughterOriginal2.GetPphi(),
       DaughterOriginal2.GetR(),DaughterOriginal2.GetTheta(),DaughterOriginal2.GetPhi());
}
            //InfoTuple->Fill(CONTAINER);
            if(RelMomD_CmDaughters<ClosePairLimit) InfoTuple_ClosePairs->Fill(CONTAINER);
        }


        KittyEvent->Reset();
        KittyFilteredEvent->Reset();

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }

    }//while(!feof(InFile))

    if(ProgressBar){
        printf("\r\033[K");
    }



    //h_rkAngle_Mom2->GetXaxis()->SetTitle("Angle (rad)");
    //h_rkAngle_Mom2->GetYaxis()->SetTitle("p_{res} (MeV)");

    //h_rkAngle_Mom2->Write();
    //h_kkAngle_Mom2->Write();

    //h_rkAngle_k_Reso->Write();
    //h_rkAngle_Reso->Write();
    //h_rkAngle_k_Daughters->Write();
    //h_rkAngle_k_Prim->Write();
    //h_rkAngle_Mom2->Write();
    //h_rkAngle_MomReso2->Write();
    //h_rkAngle_MomResoA->Write();
    //h_rkAngle_MomResoB->Write();
    h_AP_D->Write();
    h_AP_D_Large_AP_P1->Write();
    h_AP_D_Large_AP_P2->Write();
    h_AP_D_LargeAngles->Write();
    h1_AP_D->Write();
    h_AD_D->Write();
    h_AP_D1->Write();
    h_AD_D1->Write();
    h_AP_D2->Write();
    h_AD_D2->Write();
    h_AP_P1->Write();
    h_AP_P1_AR_P1P2->Write();
    h_AP_P1_dR_P1P2->Write();
    h_AP_P1_LargeAngles->Write();
    h_AD_P1->Write();
    h_AP_P2->Write();
    h_AD_P2->Write();
    h_MomDist_P1->Write();
    h_MomDist_P2->Write();
    h_MomDist_P1P2->Write();
    h_TauDist_P1->Write();
    h_TauDist_P2->Write();
    h_TauDist_P1P2->Write();
    h_AP_LAB->Write();
    h_AP_LAB_AR_LAB->Write();
    h_AP_LAB_AR_LAB_Large_AP_P1->Write();
    h_AP_LAB_dR_LAB->Write();
    h_AngleRCP1_D->Write();
    h_AngleRCP2_D->Write();
    h_AngleP1P2_D->Write();
    h_AngleRCP1_RCP2_D->Write();
    h_AngleRCP1_P1P2_D->Write();
    h_AngleRCP2_P1P2_D->Write();
    hResoMass->Write();
    hNbody->Write();
    //InfoTuple->Write();
    InfoTuple_ClosePairs->Write();

    delete KittyEvent;
    delete KittyFilteredEvent;

    delete h_AP_D;
    delete h_AP_D_Large_AP_P1;
    delete h_AP_D_Large_AP_P2;
    delete h_AP_D_LargeAngles;
    delete h1_AP_D;
    delete h_AD_D;
    delete h_AP_D1;
    delete h_AD_D1;
    delete h_AP_D2;
    delete h_AD_D2;
    delete h_AP_P1;
    delete h_AP_P1_LargeAngles;
    delete h_AD_P1;
    delete h_AP_P2;
    delete h_AD_P2;

    //delete h_rkAngle_k_Daughters;
    //delete h_rkAngle_k_Prim;
    //delete h_rkAngle_MomReso2;
    //delete h_rkAngle_MomResoA;
    //delete h_rkAngle_MomResoB;
    //delete h_rkAngle_Mom2;
    delete h_MomDist_P1;
    delete h_MomDist_P2;
    delete h_MomDist_P1P2;
    delete h_TauDist_P1;
    delete h_TauDist_P2;
    delete h_TauDist_P1P2;
    delete h_AP_LAB;
    delete h_AP_P1_AR_P1P2;
    delete h_AP_P1_dR_P1P2;
    delete h_AP_LAB_AR_LAB;
    delete h_AP_LAB_AR_LAB_Large_AP_P1;
    delete h_AP_LAB_dR_LAB;
    delete h_AngleRCP1_D;
    delete h_AngleRCP2_D;
    delete h_AngleP1P2_D;
    delete h_AngleRCP1_RCP2_D;
    delete h_AngleRCP1_P1P2_D;
    delete h_AngleRCP2_P1P2_D;
    delete hResoMass;
    delete hNbody;
    delete InfoTuple;
    delete InfoTuple_ClosePairs;
    //delete h_kkAngle_Mom2;
    //delete eventRan;
    delete eventRan1;
    delete eventRan2;
    //delete h_rkAngle_Reso;
    //delete h_rkAngle_k_Reso;
    delete fOut;
}

//we produce some particle at an angle theta
//they fly
void TestVectorOrientation(){

}

//Descr2 system
//Descr3 Core Corona All
void CorrFromEpos(TString Descr2, TString Descr3, const double& CoronaFrac){

    const unsigned NumBlankHeaderLines=3;
    const unsigned MaxPairsToRead = 16e6;
    const unsigned HighMultLimit = 128;
    const double FractionForME = 0.002;

    unsigned TotNumEvents=0;
    unsigned RejectedHighMultEvents=0;

    //bool ParticleAntiParticle = true;

    double minMom[2];
    minMom[0] = 0.4;
    minMom[1] = 0.4;

    double etaCut[2];
    etaCut[0] = 0.8;
    etaCut[1] = 0.8;

    int pdgID[2];
    pdgID[0]=2212;
    pdgID[1]=-2212;

    int eposID[2];
    eposID[0]=1120;
    eposID[1]=-1120;

    if(Descr2.Contains("pp")){
        pdgID[0]=2212;
        eposID[0]=1120;
        minMom[0] = 0.4;
        etaCut[0] = 0.8;

        pdgID[1]=2212;
        eposID[1]=1120;
        minMom[1] = 0.4;
        etaCut[1] = 0.8;
    }
    if(Descr2.Contains("pL")){
        pdgID[0]=2212;
        eposID[0]=1120;
        minMom[0] = 0.4;
        etaCut[0] = 0.8;

        pdgID[1]=3122;
        eposID[1]=2130;
        minMom[1] = 0.3;
        etaCut[1] = 0.8;
    }
    if(Descr2.Contains("LL")){
        pdgID[0]=3122;
        eposID[0]=2130;
        minMom[0] = 0.4;
        etaCut[0] = 0.8;

        pdgID[1]=3122;
        eposID[1]=2130;
        minMom[1] = 0.3;
        etaCut[1] = 0.8;
    }
    if(Descr2.Contains("pAp")){
        pdgID[0]=2212;
        eposID[0]=1120;
        minMom[0] = 0.4;
        etaCut[0] = 0.8;

        pdgID[1]=-2212;
        eposID[1]=-1120;
        minMom[1] = 0.4;
        etaCut[1] = 0.8;
    }
    if(Descr2.Contains("pAL")){
        pdgID[0]=2212;
        eposID[0]=1120;
        minMom[0] = 0.4;
        etaCut[0] = 0.8;

        pdgID[1]=-3122;
        eposID[1]=-2130;
        minMom[1] = 0.3;
        etaCut[1] = 0.8;
    }
    if(Descr2.Contains("LAL")){
        pdgID[0]=3122;
        eposID[0]=2130;
        minMom[0] = 0.4;
        etaCut[0] = 0.8;

        pdgID[1]=-3122;
        eposID[1]=-2130;
        minMom[1] = 0.3;
        etaCut[1] = 0.8;
    }

    if(Descr3.Contains("Corona")){
        pdgID[0]=eposID[0];
        pdgID[1]=eposID[1];
    }
    if(Descr3.Contains("Core")){
        eposID[0]=pdgID[0];
        eposID[1]=pdgID[1];
    }

    const TString FileToCompareWith = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/ALL_CorrectedMC_EXP/fOutput_MaybeImportant.root";
    const TString HistoToCompareWith = "hOriginal_pAp";
    TFile* fCompare = new TFile(FileToCompareWith, "read");
    TH1F* hcOriginal = (TH1F*)fCompare->Get(HistoToCompareWith);
    DLM_Histo<float> dlmOriginal;
    dlmOriginal.SetUp(1);
    dlmOriginal.SetUp(0,hcOriginal->GetNbinsX(),hcOriginal->GetBinLowEdge(1),hcOriginal->GetXaxis()->GetBinUpEdge(hcOriginal->GetNbinsX()));
    dlmOriginal.Initialize();
//printf("NB = %i %f %f\n",hcOriginal->GetNbinsX(),hcOriginal->GetBinLowEdge(1),hcOriginal->GetXaxis()->GetBinUpEdge(hcOriginal->GetNbinsX()));

    for(unsigned uBin=0; uBin<dlmOriginal.GetNbins(); uBin++){
        dlmOriginal.SetBinContent(uBin,hcOriginal->GetBinContent(uBin+1));
//printf(" %u = %f\n",);
        dlmOriginal.SetBinError(uBin,hcOriginal->GetBinError(uBin+1));
    }
    delete fCompare;

    TH1F* hCompareCk = new TH1F("hCompareCk","hCompareCk",256,0,8192);
    TH1F* hcEpos_Ck;// = new TH1F("hcEpos_Ck","hcEpos_Ck",256,0,8192);
    TH1F* hcEpos_CdPhi;// = new TH1F("hcEpos_CdPhi","hcEpos_CdPhi",256,-3.15*2,3.15*2);

    TH1F* hEpos_CkSe = new TH1F("hEpos_CkSe","hEpos_CkSe",256,0,8192);
    TH1F* hEpos_CkMe = new TH1F("hEpos_CkMe","hEpos_CkMe",256,0,8192);
    TH1F* hEpos_Ck = new TH1F("hEpos_Ck","hEpos_Ck",256,0,8192);
    TH1F* hEpos_Ck_SillyNorm;

    TH1F* hEpos_CdPhiSe = new TH1F("hEpos_CdPhiSe","hEpos_CdPhiSe",256,-3.15*2,3.15*2);
    TH1F* hEpos_CdPhiMe = new TH1F("hEpos_CdPhiMe","hEpos_CdPhiMe",256,-3.15*2,3.15*2);
    TH1F* hEpos_CdPhi = new TH1F("hEpos_CdPhi","hEpos_CdPhi",256,-3.15*2,3.15*2);

    for(unsigned uBin=1; uBin<=hCompareCk->GetNbinsX(); uBin++){
        double Momentum = hCompareCk->GetBinCenter(uBin);
        hCompareCk->SetBinContent(uBin,dlmOriginal.Eval(&Momentum));
        hCompareCk->SetBinError(uBin,dlmOriginal.EvalError(&Momentum));
    }

    CatsParticle KittyParticle;
    //CatsParticle KittyParticleCM;

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;
    TRandom3 rangen(11);

    TString InputFileName;
    //if(Descr2=="pp_All") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_All.f19";
    //if(Descr2=="All") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_All.f19";
    //if(Descr2=="Core") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_Core.f19";
    //if(Descr2=="Corona") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_Corona.f19";
    //if(Descr2=="pAp_All") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAL_All.f19";
    //if(Descr2=="pAp_Core") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAL_All.f19";
    //if(Descr2=="pAp_Corona") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAL_All.f19";
    //if(Descr2=="pAL_All") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAL_All.f19";
    //if(Descr2=="pAL_Core") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAL_All.f19";
    //if(Descr2=="pAL_Corona") InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAL_All.f19";
    InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAL_All.f19";
    //InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/ForVale/EPOS_20200121/pAp_Corona.f19";

    unsigned NumTotalPairs=0;

    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    double Time;
    int pTotal;
    int pTotalOld;
    float ProgressLoad;
    DLM_Timer dlmTimer;
    bool ProgressBar=false;

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
        return;
    }
    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;
    char* cdummy = new char [512];
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 511, InFile)){
            printf("Issue!\n");
            continue;
        }
    }
    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    CatsEvent KittyMixedEvent(pdgID[0],pdgID[1]);
    CatsEvent KittyMixedEventAp(-pdgID[0],-pdgID[1]);

    CatsEvent KittyMixedEventCM(pdgID[0],pdgID[1]);
    CatsEvent KittyMixedEventApCM(-pdgID[0],-pdgID[1]);

    while(!feof(InFile)){
//printf("Next event...\n");
        if(NumTotalPairs>=MaxPairsToRead) break;
        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("Some fscanf issue!\n");
            continue;
        }
//printf("Go on...\n");
        TotNumEvents++;
        if(NumPartInEvent>int(HighMultLimit)) RejectedHighMultEvents++;

        CatsEvent KittyEvent(pdgID[0],pdgID[1]);
        CatsEvent KittyEventAp(-pdgID[0],-pdgID[1]);
        CatsEvent KittyEventCM(pdgID[0],pdgID[1]);
        CatsEvent KittyEventApCM(-pdgID[0],-pdgID[1]);
//printf("Go to the loop...\n");
        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            //KittyParticleCM = KittyParticle;

            if(NumTotalPairs>=MaxPairsToRead) continue;
            if(NumPartInEvent>int(HighMultLimit)) continue;

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("WARNING! Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }
//printf("Original pid = %i\n",KittyParticle.GetPid());
            //reduce the corona fraction
            if(KittyParticle.GetPid()==eposID[0]){
                if(Descr2.Contains("All")&&rangen.Uniform()>CoronaFrac) continue;
                KittyParticle.SetPid(pdgID[0]);
            }
            if(KittyParticle.GetPid()==eposID[1]){
                if(Descr2.Contains("All")&&rangen.Uniform()>CoronaFrac) continue;
                KittyParticle.SetPid(pdgID[1]);
            }
//printf("Updated pid = %i\n",KittyParticle.GetPid());
            bool pid1_compatible = (KittyEvent.GetPidParticle1()==KittyParticle.GetPid()||KittyEventAp.GetPidParticle1()==KittyParticle.GetPid());
            bool pid2_compatible = (KittyEvent.GetPidParticle2()==KittyParticle.GetPid()||KittyEventAp.GetPidParticle2()==KittyParticle.GetPid());
//printf(" pid1_compatible = %i\n",pid1_compatible);
//printf(" pid2_compatible = %i\n",pid2_compatible);
//printf(" KittyEvent.GetPidParticle1() = %i\n",KittyEvent.GetPidParticle1());
//printf(" KittyEvent.GetPidParticle2() = %i\n",KittyEvent.GetPidParticle2());
//printf(" KittyEventAp.GetPidParticle1() = %i\n",KittyEventAp.GetPidParticle1());
//printf(" KittyEventAp.GetPidParticle2() = %i\n",KittyEventAp.GetPidParticle2());
            if(!pid1_compatible&&!pid2_compatible) continue;

            ///ALICE ACCEPTANCE
            if(pid1_compatible&&KittyParticle.GetP()<minMom[0]) continue;
            if(pid2_compatible&&KittyParticle.GetP()<minMom[1]) continue;
            if(pid1_compatible&&fabs(KittyParticle.GetPseudoRap())>etaCut[0]) continue;
            if(pid2_compatible&&fabs(KittyParticle.GetPseudoRap())>etaCut[1]) continue;

//printf(" pid = %i\n",KittyParticle.GetPid());
            bool TakeForMixing = FractionForME>rangen.Uniform();

            if(KittyParticle.GetPid()==KittyEvent.GetPidParticle1() || KittyParticle.GetPid()==KittyEvent.GetPidParticle2()){
                KittyEvent.AddParticle(KittyParticle);
                KittyEventCM.AddParticle(KittyParticle);
                if(TakeForMixing){
                    KittyMixedEvent.AddParticle(KittyParticle);
                    KittyMixedEventCM.AddParticle(KittyParticle);
                }
            }
            else if( (KittyParticle.GetPid()==KittyEventAp.GetPidParticle1() || KittyParticle.GetPid()==KittyEventAp.GetPidParticle2()) &&
                      KittyEventAp.GetPidParticle1()!=-KittyEventAp.GetPidParticle2() ){
                KittyEventAp.AddParticle(KittyParticle);
                KittyEventApCM.AddParticle(KittyParticle);
                if(TakeForMixing){
                    KittyMixedEventAp.AddParticle(KittyParticle);
                    KittyMixedEventApCM.AddParticle(KittyParticle);
                }
            }
            else{
                printf("Unexpected continue call\n");
                continue;
            }

        }
//printf(" Computing...\n");
        KittyEvent.ComputeParticlePairs(false,false);
        KittyEventAp.ComputeParticlePairs(false,false);
        KittyEventCM.ComputeParticlePairs(false,true);
        KittyEventApCM.ComputeParticlePairs(false,true);
//printf(" Computed...\n");
        for(unsigned uPart=0; uPart<KittyEvent.GetNumPairs(); uPart++){
//if(uPart==0) printf("PP = %u\n",KittyEvent.GetNumPairs());
            hEpos_CkSe->Fill(KittyEventCM.GetParticlePair(uPart).GetP()*500.);
            hEpos_Ck->Fill(KittyEventCM.GetParticlePair(uPart).GetP()*500.);
            hEpos_CdPhiSe->Fill(KittyEvent.GetParticlePair(uPart).GetParticle(0).GetPhi()-KittyEvent.GetParticlePair(uPart).GetParticle(1).GetPhi());
            hEpos_CdPhi->Fill(KittyEvent.GetParticlePair(uPart).GetParticle(0).GetPhi()-KittyEvent.GetParticlePair(uPart).GetParticle(1).GetPhi());
        }
        for(unsigned uPart=0; uPart<KittyEventAp.GetNumPairs(); uPart++){
//if(uPart==0) printf("PAP = %u\n",KittyEventAp.GetNumPairs());
            hEpos_CkSe->Fill(KittyEventApCM.GetParticlePair(uPart).GetP()*500.);
            hEpos_Ck->Fill(KittyEventApCM.GetParticlePair(uPart).GetP()*500.);
            hEpos_CdPhiSe->Fill(KittyEventAp.GetParticlePair(uPart).GetParticle(0).GetPhi()-KittyEventAp.GetParticlePair(uPart).GetParticle(1).GetPhi());
            hEpos_CdPhi->Fill(KittyEventAp.GetParticlePair(uPart).GetParticle(0).GetPhi()-KittyEventAp.GetParticlePair(uPart).GetParticle(1).GetPhi());
        }

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }
//printf("Next please...\n");
    }//while(!feof(InFile))
//printf("Mixed event...\n");

    KittyMixedEvent.ComputeParticlePairs(false,false);
    KittyMixedEventAp.ComputeParticlePairs(false,false);
    KittyMixedEventCM.ComputeParticlePairs(false,true);
    KittyMixedEventApCM.ComputeParticlePairs(false,true);

    for(unsigned uPart=0; uPart<KittyMixedEvent.GetNumPairs(); uPart++){
        hEpos_CkMe->Fill(KittyMixedEventCM.GetParticlePair(uPart).GetP()*500.);
        hEpos_CdPhiMe->Fill(KittyMixedEvent.GetParticlePair(uPart).GetParticle(0).GetPhi()-KittyMixedEvent.GetParticlePair(uPart).GetParticle(1).GetPhi());
    }
    for(unsigned uPart=0; uPart<KittyMixedEventAp.GetNumPairs(); uPart++){
        hEpos_CkMe->Fill(KittyMixedEventApCM.GetParticlePair(uPart).GetP()*500.);
        hEpos_CdPhiMe->Fill(KittyMixedEventAp.GetParticlePair(uPart).GetParticle(0).GetPhi()-KittyMixedEventAp.GetParticlePair(uPart).GetParticle(1).GetPhi());
    }

    hEpos_Ck_SillyNorm = (TH1F*)hEpos_CkSe->Clone("hEpos_Ck_SillyNorm");
    hEpos_CkMe->Scale(1./hEpos_CkMe->Integral(hEpos_CkMe->FindBin(180.),hEpos_CkMe->FindBin(280.)),"width");
    hEpos_Ck_SillyNorm->Scale(1./hEpos_Ck_SillyNorm->Integral(hEpos_Ck_SillyNorm->FindBin(180.),hEpos_Ck_SillyNorm->FindBin(280.)),"width");
    hEpos_Ck_SillyNorm->Divide(hEpos_CkMe);

    hEpos_CkSe->Scale(1./hEpos_CkSe->Integral(),"width");
    hEpos_CkMe->Scale(1./hEpos_CkMe->Integral(),"width");
    hEpos_Ck->Scale(1./hEpos_Ck->Integral(),"width");

    hEpos_CdPhiSe->Scale(1./hEpos_CdPhiSe->Integral(),"width");
    hEpos_CdPhiMe->Scale(1./hEpos_CdPhiMe->Integral(),"width");
    hEpos_CdPhi->Scale(1./hEpos_CdPhi->Integral(),"width");

    hEpos_Ck->Divide(hEpos_CkMe);
    hEpos_CdPhi->Divide(hEpos_CdPhiMe);

    hcEpos_Ck = (TH1F*)hEpos_Ck->Clone("hcEpos_Ck");
    //hcEpos_CdPhi = (TH1F*)hEpos_CdPhi->Clone("hcEpos_CdPhi");

    hcEpos_Ck->Divide(hCompareCk);


    TFile* fOutput;
    if(!Descr2.Contains("All")){
        fOutput = new TFile(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/CorrFromEpos/fOutput_%s_%s.root",Descr2.Data(),Descr3.Data()),"recreate");
    }
    else{
        fOutput = new TFile(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/CorrFromEpos/fOutput_%s_%s_%.2f.root",Descr2.Data(),Descr3.Data(),CoronaFrac),"recreate");
    }

    hEpos_CkSe->Write();
    hEpos_CkMe->Write();
    hEpos_Ck->Write();
    hEpos_Ck_SillyNorm->Write();
    hCompareCk->Write();
    hcEpos_Ck->Write();

    hEpos_CdPhiSe->Write();
    hEpos_CdPhiMe->Write();
    hEpos_CdPhi->Write();



    delete hEpos_CkSe;
    delete hEpos_CkMe;
    delete hEpos_Ck;
    delete hcEpos_Ck;
    delete hCompareCk;
    delete hEpos_CdPhiSe;
    delete hEpos_CdPhiMe;
    delete hEpos_CdPhi;

    delete fOutput;

}


int MIXEDEVENTS(int argc, char *argv[]){
    //DifferentTechniquesTest1("DimiPhi","pp");
    //DifferentTechniquesTest1("DimiPhi","pLambda");

    //DifferentTechniquesTest1("ForVale","pAp","Core");
    //DifferentTechniquesTest1("Bhawani","pp","with");
    DifferentTechniquesTest1("Bhawani","pp","without");
return 0;

    //ReferenceSampleStudy_2("ForMax","pi_pi");
    //ReferenceSampleStudy_2("ForMax","pi_piReso");
    //ReferenceSampleStudy_2("ForMax","piReso_piReso");
    //ReferenceSampleStudy_2("ForMax","p_KaonReso");
    //ReferenceSampleStudy_2("ForMax","pReso_KaonReso");
    //ReferenceSampleStudy_2("ForMaxRamona","pi_SigReso");
    //ReferenceSampleStudy_2("ForMaxRamona","piReso_Sig");
    //ReferenceSampleStudy_2("ForMaxRamona","piReso_SigReso");
    if(atoi(argv[1])==1) ReferenceSampleStudy_2("ForRamona","p_KaonReso");
    else if(atoi(argv[1])==2) ReferenceSampleStudy_2("ForRamona","pReso_Kaon");
    else ReferenceSampleStudy_2("ForRamona","pReso_KaonReso");
    return 0;


    //SE_for_Vale("pL");


    //CorrFromEpos("Core",1);
    //CorrFromEpos("Corona",1);
    //CorrFromEpos("All",1.0);
    //CorrFromEpos("pp","Core",1.00);
    //CorrFromEpos("pp","Corona",1.00);
    //CorrFromEpos("pp","All",1.00);
    //CorrFromEpos("pAp","Core",1.00);
    //CorrFromEpos("pAp","Corona",1.00);
    //CorrFromEpos("pAp","All",1.00);
    //CorrFromEpos("pL","Core",1.00);
    //CorrFromEpos("pL","Corona",1.00);
    //CorrFromEpos("pL","All",1.00);
    //CorrFromEpos("pAL","Core",1.00);
    //CorrFromEpos("pAL","Corona",1.00);
    //CorrFromEpos("pAL","All",1.00);
    //CorrFromEpos("LL","Core",1.00);
    //CorrFromEpos("LL","Corona",1.00);
    //CorrFromEpos("LL","All",1.00);
    //CorrFromEpos("LAL","Core",1.00);
    //CorrFromEpos("LAL","Corona",1.00);
    //CorrFromEpos("LAL","All",1.00);

    //CorrFromEpos("All",0.8);
    //CorrFromEpos("All",0.7);
    //CorrFromEpos("All",0.6);
    //CorrFromEpos("All",0.5);
    //CorrFromEpos("All",0.4);
    //CorrFromEpos("All",0.39);
    //CorrFromEpos("All",0.38);
    //CorrFromEpos("All",0.37);
    //CorrFromEpos("All",0.36);
    //CorrFromEpos("All",0.35);
    //CorrFromEpos("All",0.34);
    //CorrFromEpos("All",0.33);
    //CorrFromEpos("All",0.32);
    //CorrFromEpos("All",0.31);
    //CorrFromEpos("All",0.3);
    //CorrFromEpos("All",0.28);
    //CorrFromEpos("All",0.26);
    //CorrFromEpos("All",0.24);
    //CorrFromEpos("All",0.22);
    //CorrFromEpos("All",0.2);
    //CorrFromEpos("All",0.18);
    //CorrFromEpos("All",0.16);
    //CorrFromEpos("All",0.14);
    //CorrFromEpos("All",0.12);
    //CorrFromEpos("All",0.1);

    //ReferenceSampleStudy_1("DimiPhi","pOmega");
    //ReferenceSampleStudy_1("DimiPhi","pReso");
    //ReferenceSampleStudy_2("DimiPhi","p_p");
    //ReferenceSampleStudy_2("EposDisto","p_pReso");
    //ReferenceSampleStudy_2("EposDisto","pReso_pReso");
    //ReferenceSampleStudy_2("EposDisto","p_LamReso");
    //ReferenceSampleStudy_2("EposDisto","pReso_Lam");
    //ReferenceSampleStudy_2("EposDisto","pReso_LamReso");
    //ReferenceSampleStudy_2("EposDisto","pReso_Xim");
    //ReferenceSampleStudy_2("EposDisto","pReso_Omega");
    //ReferenceSampleStudy_2("DimiPhi","p_p");
    //ReferenceSampleStudy_2("EposDisto","Lam_Lam");
    //ReferenceSampleStudy_2("EposDisto","Lam_LamReso");
    //ReferenceSampleStudy_2("EposDisto","LamReso_LamReso");

    //ReferenceSampleStudy_2("Epos2body","p_pReso");
    //ReferenceSampleStudy_2("Only3body","p_pReso_3body");
    //ReferenceSampleStudy_2("Only3body","pReso_pReso_3body");
    //ReferenceSampleStudy_2("Epos","pReso_Phi");
    //ReferenceSampleStudy_2("Epos","LamReso_Xim");

    dEta_dPhi_Ck_QS("QS", "pp", true);
    //CompareReferenceSamples("pp");
    //CompareSameMixedEventToBoltzmann();
    //RatioBetweenLevy();
    //Fit_pL_Splines();
    //updated version of the previous one, done from scratch to avoid bugs
    //DifferentTechniquesVer1("pp");
}
