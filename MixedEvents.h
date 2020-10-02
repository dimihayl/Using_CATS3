
#ifndef MIXEDEVENTS_H
#define MIXEDEVENTS_H

#include "TString.h"

const TString TransportFile_pp_Hades =
"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CppLocal/UrQMD_Crab_1/InputFiles/"
"BigProtonSample/BigProtonSample_FM3.f19";
const TString TransportFile_pL_Hades =
"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CppLocal/UrQMD_Crab_1/InputFiles/"
"BigLambdaSample/BigLambdaSample_pL_full.f19";
const TString TransportFile_pp_Alice_small =
"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_14641LrzProtonFiles.f19";
const TString TransportFile_pp_Alice =
"/home/dimihayl/Mount/nx3/scratch9/dmihaylov/Temp/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19";
const TString TransportFile_pp_AliceCore100Corona50 =
"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_16KiLrzProtonFiles_Core100Corona50.f19";
const TString TransportFile_pp_AliceCore0Corona100 =
"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_8KiLrzProtonFiles_Core0Corona100.f19";
const TString TransportFile_pp_AliceCore100Corona0 =
"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_8KiLrzProtonFiles_Core100Corona0.f19";
const TString TransportFile_pL_Alice =
"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19";


int MIXEDEVENTS(int narg, char** ARGS);

class CatsParticle;
class CatsEvent;
class CatsDataBuffer;


class ReferenceSample{
public:
private:

};


#endif // MIXEDEVENTS_H
