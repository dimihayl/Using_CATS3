#ifndef FUNWITHCECA_H
#define FUNWITHCECA_H

class TString;

void CECA_sim_for_AI_pp_v1(int SEED, unsigned NUM_CPU, unsigned NUM_CONFIGS, unsigned K_YIELD_PER_CONFIG, int flag, TString OutputFolderName);
void CECA_sim_for_AI_pp_v2(int SEED, unsigned NUM_CPU, unsigned NUM_CONFIGS, unsigned K_YIELD_PER_CONFIG, 
    float min_rd, float max_rd, float min_ht, float max_ht, float min_tau, float max_tau,
    int flag, TString OutputFolderName);

int FUN_WITH_CECA(int narg, char** ARGS);

#endif
