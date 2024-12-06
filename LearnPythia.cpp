#include "Pythia8/Pythia.h" // access to Pythia objects.
using namespace Pythia8; // allow simplified notation.

#include <unistd.h>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

#include "DLM_Histo.h"
#include "CATStools.h"
#include "DLM_RootWrapper.h"
#include "EnvVars.h"

/*
bool CommonGrandGrandParent(Pythia8::Event& event, int id1, int id2){

}

bool CommonGrandParent(Pythia8::Event& event, int id1, int id2){

}

bool CommonParent(Pythia8::Event& event, int id1, int id2){

}
*/

struct PairPid {
    int first;
    int second;
    std::string name;
    float fraction=1;
};

//what I would like to do:
//Group the protons into pairs
//Check if each two protons have a common Ancestor
//Check if each two protons have a common Mother
//plot the k* disto for the 2 cases above, as well as for the case of no common things at all
void online_example_1(int seed, char* OutputFolder){
    // --- Initialization ---
    Pythia pythia; // Define Pythia object.
    Pythia8::Event& event = pythia.event; // quick access to current event.
    // Read in settings
    //pythia.readString("..."); // line by line...
    //pythia.readFile("cardfile.cmnd"); // or via file.
    //pythia.readString("Top:gg2ttbar = on"); // Switch on process., originally I had it on
    pythia.readString("Beams:eCM = 13000."); // 13 TeV CM energy.
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Random:setSeed = on");
    // needed for production vertex
    pythia.readString("Fragmentation:setVertices = on");
    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString(Form("Tune:pp = 14"));//Monash 
    pythia.readString("SoftQCD:all = on");    // Switch on process., originally I had it off

    pythia.readString(TString::Format("Random:seed %d", seed).Data());

    TRandom3 rangen(seed);

    int num_ca=0;
    int num_cm=0;
    int num_cg=0;

    //Pythia8::Particle p1;

    // Define histograms, external links,
    // local variables etc. here. E.g.
    //per core, c.a. 100 events / sec // The number of events to run. 1/2 hour is around 120K
    //NEW est for pp pL and pipi, we need 1.5 mins for 5K events, => 200K / hour

//100-124 :     20K x 25 = 0.5M
//125-149:     200K x 25 = 5M
//150-174:     200K x 25 = 5M
//175-200:     800K x 25 = 20M
//200-224:     800K x 25 = 20M
//225-249:       2M x 25 = 50M
    int maxEvents = 1000*100;
   // maxEvents = 10;

    const int NumMtBins = 192;
    const double MinMt = 130;
    const double MaxMt = 4970;
    //double* RangeMt = new double[NumMtBins+1];
    //RangeMt[0] = 940;
    //RangeMt[1] = 1200;
    //RangeMt[2] = 1450;
    //RangeMt[3] = 1900;
    //RangeMt[4] = 4000;

    const int NumKstarBins = 768/2;
    const float kMin = 0;
    const float kMax = 9216/2;
    const int NumRadBins = 2048;
    const float rMin = 0;
    const float rMax = 16;
    std::vector<PairPid> particle_pairs_pid;

    std::vector<int> particles_pid;
    particles_pid.push_back(2212);
    particles_pid.push_back(3122);
    particles_pid.push_back(211);

    PairPid pair_pp;
    pair_pp.first = 2212;
    pair_pp.second = 2212;
    pair_pp.name = "pp";
    particle_pairs_pid.push_back(pair_pp);

    PairPid pair_pL;
    pair_pL.first = 2212;
    pair_pL.second = 3122;
    pair_pL.name = "pL";
    particle_pairs_pid.push_back(pair_pL);

    PairPid pair_pipi;
    pair_pipi.first = 211;
    pair_pipi.second = 211;
    pair_pipi.name = "pipi";
    pair_pipi.fraction = 0.01;//originally: 0.02
    particle_pairs_pid.push_back(pair_pipi);

    const int NumPairTypes = particle_pairs_pid.size();

    const float EtaCut = 0.8;
    const float PtMin = 400;//its 300 for tha lambdas, and 500 for the protons in the analysis note (7 TeV)
    const int EventSizeMin = 0;
    const int EventSizeMax = 1e6;
    const float MtMin = 1000;
    const float MtMax = 1e6;
    const float KstarCut = 1e6;


    Int_t ComRel;
    Int_t Type;
    Int_t pid1;
    Int_t pid2;
    Float_t kstar;
    Float_t rstar;
    Float_t pt1;
    Float_t pt2;
    Float_t mT;
    TTree* pythiaTree = new TTree("pythiaTree","pythiaTree");
    pythiaTree->Branch("pid1", &pid1, "pid1/I");
    pythiaTree->Branch("pid2", &pid2, "pid2/I");    
    pythiaTree->Branch("ComRel", &ComRel, "ComRel/I");//
    pythiaTree->Branch("Type", &Type, "Type/I");//
    pythiaTree->Branch("kstar", &kstar, "kstar/F");
    pythiaTree->Branch("rstar", &rstar, "rstar/F");
    pythiaTree->Branch("pt1", &pt1, "pt1/F");
    pythiaTree->Branch("pt2", &pt2, "pt2/F");
    pythiaTree->Branch("mT", &mT, "mT/F");

    
//an ntuple for the pairs:
//k*, rstar, pt1, pt2, mT, common relative, classification (Type)
//common relative = an integer, where last digit if they have a CA at all, second to last is if they share a mother, third to last is if they share a grandmother
//classification: an integer of two digits, where each refers to one of the particles. if 0, it means the particle is primordial, if 1 it comes from a resonance


    pythia.init(); // Initialize
    printf("INITIALIZED!\n");
    printf("Num of events: %i\n", maxEvents);

    DLM_Histo<float>* KstarVsMt = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsMt_CA = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsMt_CM = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsMt_CG = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsMt_NCA = new DLM_Histo<float>[NumPairTypes];

    DLM_Histo<float>* KstarVsRstar = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsRstar_CA = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsRstar_CM = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsRstar_CG = new DLM_Histo<float>[NumPairTypes];
    DLM_Histo<float>* KstarVsRstar_NCA = new DLM_Histo<float>[NumPairTypes];

    DLM_Histo<float>* SP_x = new DLM_Histo<float>[particles_pid.size()];
    DLM_Histo<float>* SP_y = new DLM_Histo<float>[particles_pid.size()];
    DLM_Histo<float>* SP_z = new DLM_Histo<float>[particles_pid.size()];
    DLM_Histo<float>* SP_time = new DLM_Histo<float>[particles_pid.size()];
    DLM_Histo<float>* SP_r = new DLM_Histo<float>[particles_pid.size()];

    for(int ip=0; ip<NumPairTypes; ip++){
        KstarVsMt[ip].SetUp(2);
        KstarVsMt[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsMt[ip].SetUp(1, NumMtBins, MinMt, MaxMt);
        KstarVsMt[ip].Initialize();
        
        KstarVsMt_CA[ip].SetUp(2);
        KstarVsMt_CA[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsMt_CA[ip].SetUp(1, NumMtBins, MinMt, MaxMt);
        KstarVsMt_CA[ip].Initialize();

        KstarVsMt_CM[ip].SetUp(2);
        KstarVsMt_CM[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsMt_CM[ip].SetUp(1, NumMtBins, MinMt, MaxMt);
        KstarVsMt_CM[ip].Initialize();

        KstarVsMt_CG[ip].SetUp(2);
        KstarVsMt_CG[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsMt_CG[ip].SetUp(1, NumMtBins, MinMt, MaxMt);
        KstarVsMt_CG[ip].Initialize();

        KstarVsMt_NCA[ip].SetUp(2);
        KstarVsMt_NCA[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsMt_NCA[ip].SetUp(1, NumMtBins, MinMt, MaxMt);
        KstarVsMt_NCA[ip].Initialize();

        KstarVsRstar[ip].SetUp(2);
        KstarVsRstar[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsRstar[ip].SetUp(1, NumRadBins,rMin,rMax);
        KstarVsRstar[ip].Initialize();
        
        KstarVsRstar_CA[ip].SetUp(2);
        KstarVsRstar_CA[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsRstar_CA[ip].SetUp(1, NumRadBins,rMin,rMax);
        KstarVsRstar_CA[ip].Initialize();

        KstarVsRstar_CM[ip].SetUp(2);
        KstarVsRstar_CM[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsRstar_CM[ip].SetUp(1, NumRadBins,rMin,rMax);
        KstarVsRstar_CM[ip].Initialize();

        KstarVsRstar_CG[ip].SetUp(2);
        KstarVsRstar_CG[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsRstar_CG[ip].SetUp(1, NumRadBins,rMin,rMax);
        KstarVsRstar_CG[ip].Initialize();

        KstarVsRstar_NCA[ip].SetUp(2);
        KstarVsRstar_NCA[ip].SetUp(0, NumKstarBins, kMin, kMax);
        KstarVsRstar_NCA[ip].SetUp(1, NumRadBins,rMin,rMax);
        KstarVsRstar_NCA[ip].Initialize();       
    }

    for(int isp=0; isp<particles_pid.size(); isp++){
        SP_x[isp].SetUp(1);
        SP_x[isp].SetUp(0,NumRadBins,-rMax*0.5,rMax*0.5);
        SP_x[isp].Initialize();

        SP_y[isp].SetUp(1);
        SP_y[isp].SetUp(0,NumRadBins,-rMax*0.5,rMax*0.5);
        SP_y[isp].Initialize();

        SP_z[isp].SetUp(1);
        SP_z[isp].SetUp(0,NumRadBins,-rMax*0.5,rMax*0.5);
        SP_z[isp].Initialize();

        SP_time[isp].SetUp(1);
        SP_time[isp].SetUp(0,NumRadBins,rMin,rMax);
        SP_time[isp].Initialize();

        SP_r[isp].SetUp(1);
        SP_r[isp].SetUp(0,NumRadBins,rMin,rMax);
        SP_r[isp].Initialize();
    }


    TH1F* hEventSize = new TH1F("hEventSize","hEventSize",5120,0,5120);

//particles, like lambda, that have a large decay position, should be treated as primary
//and also, particle with large production position should be discated (feed-down)
//it seems that only particles with 0 prod position are really primaries
//indeed, even delta decays happen at 0 prod position in pythia

    //usleep(5e6);
    // --- The event loop ---
    for(int iEvent = 0; iEvent < maxEvents; iEvent++){
        // Generate next event;
        // Produce the next event, returns true on success.
        if(!pythia.next()) {
        // Any error handling goes here.
            printf("ERROR!!!\n");
        }

        if(pythia.event.size()<EventSizeMin || pythia.event.size()>EventSizeMax) continue;

        hEventSize->Fill(pythia.event.size());

        if(iEvent%1000 == 0)
            printf(" EVENT %i: size = %i\n", iEvent, pythia.event.size());

        //printf("es%i\n",event.size()-1);
        for(int iPart1=3; iPart1<pythia.event.size()-1; iPart1++){
            //printf("%i\n",iPart1);

            pid1 = pythia.event.at(iPart1).id();
            pt1 = pythia.event.at(iPart1).pT()*1000;
            bool FoundIt = false;
            for(int ip=0; ip<particle_pairs_pid.size(); ip++){
                if( fabs(pid1) == fabs(particle_pairs_pid.at(ip).first) || fabs(pid1) == fabs(particle_pairs_pid.at(ip).second)){
                    FoundIt = true;
                }
            }
            if(!FoundIt) continue;
            if(pythia.event.at(iPart1).isHadron()==false) continue;
            //if(pythia.event.at(iPart1).isFinal()==false) continue;

           // if (pythia.event.at(iPart1).tProd()!=0 && pythia.event.at(iPart1).tProd()<1e-10){
            //printf("iPart %i: pid %i\n",iPart1,pythia.event.at(iPart1).id());
            //printf(" -> mothers:");
            //vector<int> mother_list = pythia.event.at(iPart1).motherList();
            //for(const int& im : mother_list) printf(" %i", im);
            //printf("\n");
            //printf(" -> daughters:");
            //vector<int> daughter_list = pythia.event.at(iPart1).daughterList();
            //for(const int& id : daughter_list) printf(" %i", id);
            //printf("\n");   
            //printf(" fsp: %i\n", pythia.event.at(iPart1).isFinal());
            //printf(" mom: %.2e, %.2e, %.2e, %.2e\n", pythia.event.at(iPart1).e()*1000., pythia.event.at(iPart1).px()*1000.,pythia.event.at(iPart1).py()*1000.,pythia.event.at(iPart1).pz()*1000.);
            //printf(" prd: %.2e, %.2e, %.2e, %.2e\n", pythia.event.at(iPart1).tProd(), pythia.event.at(iPart1).xProd(),pythia.event.at(iPart1).yProd(),pythia.event.at(iPart1).zProd());
            //printf(" dec: %.2e, %.2e, %.2e, %.2e\n", pythia.event.at(iPart1).tDec(), pythia.event.at(iPart1).xDec(),pythia.event.at(iPart1).yDec(),pythia.event.at(iPart1).zDec());
            //printf("----------------------------\n");
            //usleep(40e3);
            //}

            CatsLorentzVector LV1;
            LV1.Set(
                pythia.event.at(iPart1).tProd()*1e12,
                pythia.event.at(iPart1).xProd()*1e12,
                pythia.event.at(iPart1).yProd()*1e12,
                pythia.event.at(iPart1).zProd()*1e12,
                pythia.event.at(iPart1).e()*1000.,
                pythia.event.at(iPart1).px()*1000.,
                pythia.event.at(iPart1).py()*1000.,
                pythia.event.at(iPart1).pz()*1000.);
            float Eta1 = LV1.GetPseudoRap();
            if(fabs(Eta1)>EtaCut) continue;
            if(fabs(LV1.GetPt())<PtMin) continue;
            if(fabs(LV1.GetR())>1000){
                //printf("LV1 = %.2e\n", LV1.GetR());
                continue;
            }

            for(int isp=0; isp<particles_pid.size(); isp++){
                if( fabs(pythia.event.at(iPart1).id())==fabs(particles_pid.at(isp)) ){
                    SP_x[isp].Fill(pythia.event.at(iPart1).xProd()*1e12);
                    SP_y[isp].Fill(pythia.event.at(iPart1).yProd()*1e12);
                    SP_z[isp].Fill(pythia.event.at(iPart1).zProd()*1e12);
                    SP_time[isp].Fill(pythia.event.at(iPart1).tProd()*1e12);
                    SP_r[isp].Fill(LV1.GetR());
                }
            }

            std::vector<int> p1_mother = pythia.event.at(iPart1).motherList();
            std::vector<int> p1_grandma;
            bool p1_primordial = true;
            for(int& mom : p1_mother){
                std::vector<int> temp_list = pythia.event.at(mom).motherList();
                p1_grandma.insert(p1_grandma.end(), temp_list.begin(), temp_list.end());
                //these are all hadrons, excluding diquarks
                if(mom>=100 && (mom/10)%10){
                    p1_primordial = false;
                }
            }
            


            //this makes it so to remove common stuff
            for(int iPart2=iPart1+1; iPart2<pythia.event.size(); iPart2++){
                //printf(" %i\n",iPart2);

                pid2 = pythia.event.at(iPart2).id();
                pt2 = pythia.event.at(iPart2).pT()*1000;

                int WhichPair = -1;
                for(int ip=0; ip<particle_pairs_pid.size(); ip++){
                    if( pid1==particle_pairs_pid.at(ip).first && pid2==particle_pairs_pid.at(ip).second ) WhichPair = ip;
                    if( pid2==particle_pairs_pid.at(ip).first && pid1==particle_pairs_pid.at(ip).second ) WhichPair = ip;

                    if( pid1==-particle_pairs_pid.at(ip).first && pid2==-particle_pairs_pid.at(ip).second ) WhichPair = ip;
                    if( pid2==-particle_pairs_pid.at(ip).first && pid1==-particle_pairs_pid.at(ip).second ) WhichPair = ip;

                }
                if(WhichPair<0) continue;
                if(pythia.event.at(iPart2).isHadron()==false) continue;
                //if(pythia.event.at(iPart2).isFinal()==false) continue;

                float rannum = rangen.Uniform();
                if(rannum>particle_pairs_pid.at(WhichPair).fraction) continue;


                //printf(" pid %i %i\n", pythia.event.at(iPart1).id(), pid2);
                //usleep(40e3);
                
                bool CommonAncestor = false;
                for(int iAnc=2; iAnc<iPart1; iAnc++){
                    //iAnc has to be smaller than iPart for this to work !!!
                    if(pythia.event.at(iPart1).isAncestor(iAnc) && pythia.event.at(iPart2).isAncestor(iAnc)){
                        CommonAncestor = true;
//printf("CA = %i (%i %i)\n", iAnc,iPart1,iPart2);
                        break;
                    }
                }
                
                CatsLorentzVector LV2;
                LV2.Set(
                pythia.event.at(iPart2).tProd()*1e12,
                pythia.event.at(iPart2).xProd()*1e12,
                pythia.event.at(iPart2).yProd()*1e12,
                pythia.event.at(iPart2).zProd()*1e12,
                pythia.event.at(iPart2).e()*1000.,
                pythia.event.at(iPart2).px()*1000.,
                pythia.event.at(iPart2).py()*1000.,
                pythia.event.at(iPart2).pz()*1000.);
                float Eta2 = LV2.GetPseudoRap();
                if(fabs(Eta2)>EtaCut) continue;
                if(fabs(LV2.GetPt())<PtMin) continue;
                if(fabs(LV2.GetR())>1000){
                    //printf("LV2 = %.2e\n", LV1.GetR());
                    continue;
                }

                bool CommonMother = false;
                bool CommonGrandMa = false;
                bool p2_primordial = true;
                std::vector<int> p2_mother = pythia.event.at(iPart2).motherList();
                std::vector<int> p2_grandma;
                for(int& mom : p2_mother){
                    std::vector<int> temp_list = pythia.event.at(mom).motherList();
                    p2_grandma.insert(p2_grandma.end(), temp_list.begin(), temp_list.end());
                    //these are all hadrons, excluding diquarks
                    if(mom>=100 && (mom/10)%10){
                        p2_primordial = false;
                    }
                }


                for(int& mom1 : p1_mother){
                    for(int& mom2 : p2_mother){
                        if(mom1==mom2){
                            CommonMother = true;
                            break;
                        }
                    }
                    if(CommonMother) break;
                }

                for(int& gm1 : p1_grandma){
                    for(int& gm2 : p2_grandma){
                        if(gm1==gm2){
                            CommonGrandMa = true;
                            break;
                        }
                    }
                    if(CommonGrandMa) break;
                }
/*
                printf("id1:%i(%i) id2:%i(%i) p1:%i p2:%i ca:%i cm:%i cg:%i\n",
                pythia.event.at(iPart1).id(), iPart1, pid2, iPart2,
                p1_primordial, p2_primordial, CommonAncestor, CommonMother, CommonGrandMa);

                if(CommonMother){
                    printf(" Sisters of (%i)\n ", iPart1);
                    for(int& sis : pythia.event.at(iPart1).sisterList()){
                        printf("%i ",sis);
                    }
                    printf("\n");
                }
*/
                num_ca += CommonAncestor;
                num_cm += CommonMother;
                num_cg += CommonGrandMa;

                //if(CommonAncestor){
                //    printf("id1:%i(%i) id2:%i(%i) p1:%i p2:%i ca:%i cm:%i cg:%i\n",
                //    pythia.event.at(iPart1).id(), iPart1, pid2, iPart2,
                //    p1_primordial, p2_primordial, CommonAncestor, CommonMother, CommonGrandMa);
                //}

                CatsLorentzVector LV_BOOST = LV1 + LV2;
                LV1.Boost(LV_BOOST);
                LV2.Boost(LV_BOOST);
                kstar = LV1.GetP();       
                mT = 0.5*LV_BOOST.GetMt();
                CatsLorentzVector Rel12 = LV1 - LV2;
                rstar = Rel12.GetR();
                Type = p1_primordial + p2_primordial*10;
                //verfied: if CommonMother => also CommonGrandMa
                ComRel = CommonAncestor + 10*CommonMother + 100*CommonGrandMa;
                //printf(" kstar q = %.0f %.0f\n", kstar*2, Rel12.GetP());

                //if(WhichPair==1){
                //    if(kstar<60){
                //        printf("%f %i %i\n",kstar, pythia.event.at(iPart1).id(), pid2);
                //    }
                //}

                //if(kstar<400 && rstar>0.1 && rstar<10)
                //    printf("%f %f\n",kstar,rstar);

                if(kstar>KstarCut){
                    continue;
                }

                KstarVsMt[WhichPair].Fill(kstar,mT);
                KstarVsRstar[WhichPair].Fill(kstar,rstar);
                //if(CommonAncestor || CommonMother || CommonGrandMa){
                if(CommonAncestor){
                    KstarVsMt_CA[WhichPair].Fill(kstar,mT);
                    KstarVsRstar_CA[WhichPair].Fill(kstar,rstar);
                }
                if(CommonMother){
                    KstarVsMt_CM[WhichPair].Fill(kstar,mT);
                    KstarVsRstar_CM[WhichPair].Fill(kstar,rstar);
                }
                if(CommonGrandMa){
                    KstarVsMt_CG[WhichPair].Fill(kstar,mT);
                    KstarVsRstar_CG[WhichPair].Fill(kstar,rstar);
                }
                if(!CommonAncestor && !CommonMother && !CommonGrandMa){
                    KstarVsMt_NCA[WhichPair].Fill(kstar,mT);
                    KstarVsRstar_NCA[WhichPair].Fill(kstar,rstar);
                }
                pythiaTree->Fill();
            }
        }

        /*
        for(int iPart=0; iPart<pythia.event.size(); iPart++){
            std::vector<int> mothers = pythia.event[iPart].motherList();
            std::vector<int> daughters = pythia.event[iPart].daughterList();
            if(pythia.event[iPart].id()==2212){
                //printf("  iPart: %i\n",iPart);
                for(int& mom : mothers){
                    //printf("   mother = %i\n", mom);
                }
                for(int& dau : daughters){
                    //printf("   daughter = %i\n", dau);
                }
                for(int aPart=0; aPart<pythia.event.size(); aPart++){
                    if (pythia.event[iPart].isAncestor(aPart)){
                        printf("   %i is isAncestor\n", aPart);
                    }
                }
            }  
        }
        */
        
        // Analyse event; fill histograms etc.
        //usleep(100e3);
    } // End event loop.
    // --- Calculate final statistics ---
    //pythia.stat();
    // Print histograms etc.

    TH2F** hKstarVsMt = new TH2F* [NumPairTypes];
    TH2F** hKstarVsMt_CA = new TH2F* [NumPairTypes];
    TH2F** hKstarVsMt_CM = new TH2F* [NumPairTypes];
    TH2F** hKstarVsMt_CG = new TH2F* [NumPairTypes];
    TH2F** hKstarVsMt_NCA = new TH2F* [NumPairTypes];

    TH2F** hKstarVsRstar = new TH2F* [NumPairTypes];
    TH2F** hKstarVsRstar_CA = new TH2F* [NumPairTypes];
    TH2F** hKstarVsRstar_CM = new TH2F* [NumPairTypes];
    TH2F** hKstarVsRstar_CG = new TH2F* [NumPairTypes];
    TH2F** hKstarVsRstar_NCA = new TH2F* [NumPairTypes];

    TH1F** hSP_x = new TH1F* [particles_pid.size()];
    TH1F** hSP_y = new TH1F* [particles_pid.size()];
    TH1F** hSP_z = new TH1F* [particles_pid.size()];
    TH1F** hSP_time = new TH1F* [particles_pid.size()];


    for(int ip=0; ip<NumPairTypes; ip++){
        std::string histo_name;

        histo_name = "hKstarVsMt_" + particle_pairs_pid.at(ip).name;
        KstarVsMt[ip].ComputeError();
        hKstarVsMt[ip] = Convert_DlmHisto_TH2F(&KstarVsMt[ip], histo_name.c_str());

        histo_name = "hKstarVsMt_CA_" + particle_pairs_pid.at(ip).name;
        KstarVsMt_CA[ip].ComputeError();
        hKstarVsMt_CA[ip] = Convert_DlmHisto_TH2F(&KstarVsMt_CA[ip], histo_name.c_str());

        histo_name = "hKstarVsMt_CM_" + particle_pairs_pid.at(ip).name;
        KstarVsMt_CM[ip].ComputeError();
        hKstarVsMt_CM[ip] = Convert_DlmHisto_TH2F(&KstarVsMt_CM[ip], histo_name.c_str());   

        histo_name = "hKstarVsMt_CG_" + particle_pairs_pid.at(ip).name;
        KstarVsMt_CG[ip].ComputeError();
        hKstarVsMt_CG[ip] = Convert_DlmHisto_TH2F(&KstarVsMt_CG[ip], histo_name.c_str());               

        histo_name = "hKstarVsMt_NCA_" + particle_pairs_pid.at(ip).name;
        KstarVsMt_NCA[ip].ComputeError();
        hKstarVsMt_NCA[ip] = Convert_DlmHisto_TH2F(&KstarVsMt_NCA[ip], histo_name.c_str());


        histo_name = "hKstarVsRstar_" + particle_pairs_pid.at(ip).name;
        KstarVsRstar[ip].ComputeError();
        hKstarVsRstar[ip] = Convert_DlmHisto_TH2F(&KstarVsRstar[ip], histo_name.c_str());

        histo_name = "hKstarVsRstar_CA_" + particle_pairs_pid.at(ip).name;
        KstarVsRstar_CA[ip].ComputeError();
        hKstarVsRstar_CA[ip] = Convert_DlmHisto_TH2F(&KstarVsRstar_CA[ip], histo_name.c_str());

        histo_name = "hKstarVsRstar_CM_" + particle_pairs_pid.at(ip).name;
        KstarVsRstar_CM[ip].ComputeError();
        hKstarVsRstar_CM[ip] = Convert_DlmHisto_TH2F(&KstarVsRstar_CM[ip], histo_name.c_str());

        histo_name = "hKstarVsRstar_CG_" + particle_pairs_pid.at(ip).name;
        KstarVsRstar_CG[ip].ComputeError();
        hKstarVsRstar_CG[ip] = Convert_DlmHisto_TH2F(&KstarVsRstar_CG[ip], histo_name.c_str());


        histo_name = "hKstarVsRstar_NCA_" + particle_pairs_pid.at(ip).name;
        KstarVsRstar_NCA[ip].ComputeError();
        hKstarVsRstar_NCA[ip] = Convert_DlmHisto_TH2F(&KstarVsRstar_NCA[ip], histo_name.c_str());


        //DLM_Histo<float> KstarVsMt_Ratio(KstarVsMt_CA);
        //KstarVsMt_Ratio /= KstarVsMt_NCA;
        //TH2F* hKstarVsMt_Ratio = Convert_DlmHisto_TH2F(&KstarVsMt_Ratio, "hKstarVsMt_Ratio");
    }

    for(int isp=0; isp<particles_pid.size(); isp++){
        std::string histo_name;

        histo_name = "hSP_x_" + std::to_string(particles_pid.at(isp));
        SP_x[isp].ComputeError();
        hSP_x[isp] = Convert_DlmHisto_TH1F(&SP_x[isp], histo_name.c_str());  

        histo_name = "hSP_y_" + std::to_string(particles_pid.at(isp));
        SP_y[isp].ComputeError();
        hSP_y[isp] = Convert_DlmHisto_TH1F(&SP_y[isp], histo_name.c_str());          

        histo_name = "hSP_z_" + std::to_string(particles_pid.at(isp));
        SP_z[isp].ComputeError();
        hSP_z[isp] = Convert_DlmHisto_TH1F(&SP_z[isp], histo_name.c_str());  

        histo_name = "hSP_time_" + std::to_string(particles_pid.at(isp));
        SP_time[isp].ComputeError();
        hSP_time[isp] = Convert_DlmHisto_TH1F(&SP_time[isp], histo_name.c_str());  

    }


    TString OutputFileName;
    if(OutputFolder){
        OutputFileName = TString::Format("%s/TomAndJerry_S%i.root",OutputFolder,seed);
    }
    else{
        OutputFileName = TString::Format("%s/LearnPythia/TomAndJerry_S%i.root",GetFemtoOutputFolder(),seed);
    }

    TFile fOutput(OutputFileName, "recreate");
    pythiaTree->Write();
    hEventSize->Write();

    for(int ip=0; ip<NumPairTypes; ip++){
        hKstarVsMt[ip]->Write();
        hKstarVsMt_CA[ip]->Write();
        hKstarVsMt_CM[ip]->Write();
        hKstarVsMt_CG[ip]->Write();
        hKstarVsMt_NCA[ip]->Write();
    }


    for(int ip=0; ip<NumPairTypes; ip++){
        hKstarVsRstar[ip]->Write();
        hKstarVsRstar_CA[ip]->Write();
        hKstarVsRstar_CM[ip]->Write();
        hKstarVsRstar_CG[ip]->Write();
        hKstarVsRstar_NCA[ip]->Write();        
    }

    for(int isp=0; isp<particles_pid.size(); isp++){
        hSP_x[isp]->Write();
        hSP_y[isp]->Write();
        hSP_z[isp]->Write();
        hSP_time[isp]->Write();
    }


    //printf("num_ca = %i, num_cm = %i, num_cg = %i\n",num_ca,num_cm,num_cg);

    return;
}





void Test_pReso_Omega_angulars(int seed, char* OutputFolder){
    // --- Initialization ---
    Pythia pythia; // Define Pythia object.
    Pythia8::Event& event = pythia.event; // quick access to current event.
    // Read in settings
    //pythia.readString("..."); // line by line...
    //pythia.readFile("cardfile.cmnd"); // or via file.
    //pythia.readString("Top:gg2ttbar = on"); // Switch on process., originally I had it on
    pythia.readString("Beams:eCM = 13600."); // 13 TeV CM energy.
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Random:setSeed = on");
    // needed for production vertex
    pythia.readString("Fragmentation:setVertices = on");
    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString(Form("Tune:pp = 14"));//Monash 
    pythia.readString("SoftQCD:all = on");    // Switch on process., originally I had it off

    pythia.readString(TString::Format("Random:seed %d", seed).Data());

    TRandom3 rangen(seed);

    int maxEvents = 10*1000;
   // maxEvents = 10;

    std::vector<PairPid> particle_pairs_pid;


    std::vector<int> primordial_pid;
    std::vector<int> resonance_pid;
    //particles_pid.push_back(2212);
    //particles_pid.push_back(2112);
    primordial_pid.push_back(4122);

    resonance_pid.push_back(2224);
    resonance_pid.push_back(2214);
    resonance_pid.push_back(2114);
    resonance_pid.push_back(1114);
    resonance_pid.push_back(12112);
    resonance_pid.push_back(1214); 
    resonance_pid.push_back(12116); 

    resonance_pid.push_back(22212); 
    resonance_pid.push_back(12212); 
    resonance_pid.push_back(2124); 
    resonance_pid.push_back(2126);
    resonance_pid.push_back(31214);
    resonance_pid.push_back(3124);
    resonance_pid.push_back(12216);
    resonance_pid.push_back(22112);
    resonance_pid.push_back(2126);
    resonance_pid.push_back(32224);
    resonance_pid.push_back(32112);
    resonance_pid.push_back(2216);
    resonance_pid.push_back(2222);
    resonance_pid.push_back(3225);

    //PairPid pair_pp;
    //pair_pp.first = 2212;
    //pair_pp.second = 2212;
    //pair_pp.name = "pp";
    //particle_pairs_pid.push_back(pair_pp);



    const float EtaCut = 0.8;
    const float PtMin = 400;//its 300 for tha lambdas, and 500 for the protons in the analysis note (7 TeV)
    const int EventSizeMin = 0;
    const int EventSizeMax = 1e6;
    const float MtMin = 1000;
    const float MtMax = 1e6;
    const float KstarCut = 1e6;


    Int_t ComRel;
    Int_t Type;
    Int_t pidLc;
    Int_t pidReso;
    Float_t kstar;
    Float_t rstar;
    Float_t ptLc;
    Float_t ptReso;
    Float_t mT;
    Float_t cos_Lc_Reso;
    TTree* pythiaTree = new TTree("pythiaTree","pythiaTree");
    pythiaTree->Branch("pidLc", &pidLc, "pidLc/I");
    pythiaTree->Branch("pidReso", &pidReso, "pidReso/I");    
    pythiaTree->Branch("kstar", &kstar, "kstar/F");//of the final pair, after decay Reso -> ppi
    pythiaTree->Branch("rstar", &rstar, "rstar/F");
    pythiaTree->Branch("cos_Lc_Reso", &cos_Lc_Reso, "cos_Lc_Reso/F");
    pythiaTree->Branch("ptLc", &ptLc, "ptLc/F");
    pythiaTree->Branch("ptReso", &ptReso, "ptReso/F");
    pythiaTree->Branch("mT", &mT, "mT/F");

    
//an ntuple for the pairs:
//k*, rstar, pt1, pt2, mT, common relative, classification (Type)
//common relative = an integer, where last digit if they have a CA at all, second to last is if they share a mother, third to last is if they share a grandmother
//classification: an integer of two digits, where each refers to one of the particles. if 0, it means the particle is primordial, if 1 it comes from a resonance


    pythia.init(); // Initialize
    printf("INITIALIZED!\n");
    printf("Num of events: %i\n", maxEvents);



    TH1F* hEventSize = new TH1F("hEventSize","hEventSize",5120,0,5120);


    // --- The event loop ---
    for(int iEvent = 0; iEvent < maxEvents; iEvent++){
        // Generate next event;
        // Produce the next event, returns true on success.
        if(!pythia.next()) {
        // Any error handling goes here.
            printf("ERROR!!!\n");
        }

        if(pythia.event.size()<EventSizeMin || pythia.event.size()>EventSizeMax) continue;

        hEventSize->Fill(pythia.event.size());

        std::vector<int> Lc_position;
        std::vector<int> Reso_position;
        if(iEvent%1000 == 0)
            printf(" EVENT %i: size = %i\n", iEvent, pythia.event.size());

        //printf("es%i\n",event.size()-1);
        for(int iPart1=3; iPart1<pythia.event.size()-1; iPart1++){
            //printf("%i\n",iPart1);

            int pid_tmp = pythia.event.at(iPart1).id();
            if(abs(pid_tmp)==abs(primordial_pid.at(0))){
                Lc_position.push_back(iPart1);
                printf("Lc pid = %i\n", pid_tmp);
            }
            for(int& reso_id : resonance_pid ){
                if(abs(pid_tmp)==abs(reso_id)){
                    Reso_position.push_back(iPart1);
                    //printf("Reso pid = %i\n", pid_tmp);
                }
            }

        }

    }

    TString OutputFileName;
    if(OutputFolder){
        OutputFileName = TString::Format("%s/TomAndJerry_S%i.root",OutputFolder,seed);
    }
    else{
        OutputFileName = TString::Format("%s/LearnPythia/Nibles_Lc_Reso_S%i.root",GetFemtoOutputFolder(),seed);
    }

    TFile fOutput(OutputFileName, "recreate");
    pythiaTree->Write();
    hEventSize->Write();

    //printf("num_ca = %i, num_cm = %i, num_cg = %i\n",num_ca,num_cm,num_cg);

    return;

}




int LEARN_PYTHIA(int argc, char *argv[]){
    printf("LEARN_PYTHIA, HELLO!!!\n\n");
    //Compare_pp_AV18_Haide_SingleGauss(atoi(argv[1]));

    //online_example_1(atoi(argv[1]), NULL);

    Test_pReso_Omega_angulars(11, NULL);


    return 0;
}