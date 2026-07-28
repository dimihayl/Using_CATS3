#include "PotentialStuff.h"

#include "EnvVars.h"



#include <iostream>
#include <stdio.h>
#include <string.h>

#include <unistd.h>

#include "TH1F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "DLM_Source.h"
#include "DLM_RootFit.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "DLM_Potentials.h"
#include "CommonAnaFunctions.h"

using namespace std;


void pot_strength_1(){

    //double ApproxPotStrength(double (*f)(double*), double* Pars, unsigned short usPW);
    double Pars1S0[10] = {0,0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0};//2 dummy stuff, cats style
    double strength_av18 = -ApproxPotStrength(fDlmPot, Pars1S0, 0) * Mass_p * 0.5 / hbarc / hbarc;

    double ParsSG[4] = {0,0, -31.69, 1.78};//2 dummy stuff, cats style
    double strength_sg = -ApproxPotStrength(SingleGauss, ParsSG, 0) * Mass_p * 0.5 / hbarc / hbarc;


    printf("strength_av18 = %f\n",strength_av18);
    printf("strength_sg = %f\n",strength_sg);

}

//FloatingGaussians
//we will use 3 gaussians, i.e. 9 parameters.
void generate_potentials_fg(TString OutputFolder, int NumIter, int SEED){
    double Pars1S0[10] = {0,0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0};//2 dummy stuff, cats style
    double strength_av18 = -ApproxPotStrength(fDlmPot, Pars1S0, 0) * Mass_p * 0.5 / hbarc / hbarc;
    
    const double min_ampl = 0.1;
    const double max_ampl = 2e5;

    const double min_pos = 0;
    const double max_pos = 2;

    const double min_rng = 0.05;
    const double max_rng = 4;

    double* ParsTriGauss = new double [12];
    ParsTriGauss[0] = 0;
    ParsTriGauss[1] = 1;
    ParsTriGauss[2] = 3;

    printf("strength_av18 = %f\n", strength_av18);

    TRandom3 rangen(SEED);
    for(int iIter=0; iIter<NumIter; iIter++){
        double ampl[3];
        double pos[3];
        double rng[3];
        double strength_tg = 0;
        //we want the two strngths to be of the same order of magnitude, more or less
        //MAKE THE CONDITION MEANINGFUL
        //while( fabs((strength_av18-strength_tg)/strength_av18) > 0.99  ){
            for(int iG=0; iG<3; iG++){
                //ln-uniform, c.a. from 0.1 MeV to over 1e5 MeV
                ampl[iG] = pow(TMath::E(), rangen.Uniform(log(min_ampl), log(max_ampl)));
                if( rangen.Uniform() > 0.5 ) ampl[iG] = -0.5;//equal prob for + or -
                //slight preference towards smaller numbers
                pos[iG] = rangen.Uniform(min_pos, max_pos);
                rng[iG] = sqrt(rangen.Uniform(min_rng, max_rng*max_rng));
                ParsTriGauss[3+iG*3] = ampl[iG];
                ParsTriGauss[3+iG*3+1] = pos[iG];
                ParsTriGauss[3+iG*3+2] = rng[iG];
            }
            strength_tg = -ApproxPotStrength(FloatingGaussians, ParsTriGauss, 0) * Mass_p * 0.5 / hbarc / hbarc;
        //}
        if(strength_tg>0.9)
            printf(" strength_tg = %f\n", strength_tg);
    }
    
    delete [] ParsTriGauss;
}


int POTENTIAL_STUFF(int argc, char *argv[]){
    printf("POTENTIAL_STUFF\n");

    //pot_strength_1();
    generate_potentials_fg("", 3000, 23);

    return 0;
}