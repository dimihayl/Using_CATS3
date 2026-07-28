#include "ThermalFistStuff.h"

#include <iostream>
#include <fstream>
#include <cmath>


//#include "BlastWaveModel.h"
//#include "ParticleDB.h"

#include "TString.h"


void cgpt_test1(){
/*
    // Load PDG database (adjust path if needed)
    thermalfist::ParticleDB pdg("data/PDG_2020.txt");

    // Define parameters
    double T = 0.100;         // temperature in GeV
    double rhoMax = 0.9;      // max transverse flow rapidity (tanh^{-1}(beta_T))
    double R = 10.0;          // system size in fm
    int profile_n = 1;        // flow profile index (1 = linear)

    // Create the model
    thermalfist::BlastWaveModel bwModel(T, rhoMax, R, profile_n);

    // Get particle info: pi+
    const thermalfist::ParticleType& pion = pdg.GetParticle(211);
    double m = pion.Mass();

    TString OutputFolderName = TString::Format("%s/ThermalFistStuff/",GetFemtoOutputFolder());

    // Calculate spectrum over pT
    std::ofstream out(OutputFolderName+"spectrum_pi_plus.txt");
    for (double pT = 0.05; pT < 2.0; pT += 0.05) {
        double value = bwModel.MtPtIntegrand(pT, m);
        out << pT << "\t" << value << "\n";
    }
    out.close();

    std::cout << "Spectrum saved to 'spectrum_pi_plus.txt'\n";
*/
}


int THERMAL_FIST_STUFF(int argc, char *argv[]){
    printf("THERMAL_FIST_STUFF\n");

    return 0;
}



