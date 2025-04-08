//
// Created by gst on 23/10/23.
//
#include <cmath>
#include <cstdio>
#include <vector>
#include <chopper-lib.h>
#include "choppers.h"

struct parameters {
  double path;
  double angle;
  double speed;
  double phase;
  double delay;
};

auto bifrost(double E_0, double L_0, double chopPulseOpening)
-> std::map<std::string, chopper_parameters>
//-> std::map<std::string, std::map<std::string, double>>
//std::map<std::string, double>
{
// Transferred parameters
    double chopPulseFrequencyOrder=14; // Number of chopper pulses pr moderator pulse. It will automatically be reduced when nesesary and a warning will be written in the promt.
    double chopBWPos=78;  // Distance from pulse shapping choppers to BW Chopper

/*************************************** Chopper Variables  *******************************************/
    double chopFrameOverlap1Pos= 8.530;    // Distance from moderator to first frame owerlap chopper
    double chopFrameOverlap2Pos= 14.973;    // Distance from moderator to second frame owerlap chopper

// If there is set a value of L_0, overwrite E_0 and calculate E_0 from L_0
    if (L_0>0){
        E_0=81.82/(L_0*L_0);
    }


/************************************************/
/*                  Chopper calculations                    */
/************************************************/

    double PulseHighFluxOffset=2.0e-4; // Time from T0 to high pulse.
    double ModPulseLengthHighF=2.86e-3; // width of high pulse

    double InstLength=162.0;
    double chopPulseDist= 4.41+0.032+2.0-0.1;  // Distance fro moderator to Pulse chapping chopper

    /******* Check if pulse shapping chopper opening is large enough for requested frequency or reduce frequency *******/
    if  (chopPulseFrequencyOrder*chopPulseOpening > 170.0/360.0/14.0) {
        chopPulseFrequencyOrder=floor(170.0/360.0/14.0/chopPulseOpening);
        printf(" \n \n Warning: Impossible combination of chopPulseFrequencyOrder and chopPulseOpening chosen, chopPulseFrequencyOrder reduced to: %f  \n", chopPulseFrequencyOrder);
    }

    auto lambda_1=1.0/(0.1106*sqrt(E_0));  /**** general chopper calculations **********/
    auto WavelengthBand = 1/(InstLength-chopPulseDist)/14.0/2.528e-4;
    auto lambda_0=lambda_1-WavelengthBand;
    auto v_0=3956.0/lambda_1;
    auto v_1=3956.0/lambda_0;

/***********  Pulse shaping chopper calculations **********/
    auto chopPulseOffset=(chopPulseDist/v_1+chopPulseDist/v_0)/2.0+ModPulseLengthHighF/2.0+PulseHighFluxOffset;
    auto chopPulsePhaseOffset=  (chopPulseOffset+ chopPulseOpening/2.0)*14.0*chopPulseFrequencyOrder*360.0-170.0/2.0;
    auto chopPulse2PhaseOffset= chopPulsePhaseOffset- 360.0*(chopPulseOpening*14.0*chopPulseFrequencyOrder)+170.0;

    if  (chopPulseFrequencyOrder == 0) {
        chopPulsePhaseOffset= 0;
        chopPulse2PhaseOffset= 0;
        printf(" \n \n Warning: Pulse shaping chopper parked! Setting the offsets to zero");
    }


/*********** Frame Overlap chopper calculations ******************/
    auto chopFrameOverlap1Open= 1.0/14.0/InstLength*(chopFrameOverlap1Pos)*1.5 ;
    auto chopFrameOverlap1Offset=(  ( (chopFrameOverlap1Pos)/v_1+(chopFrameOverlap1Pos)/v_0)/2.0+PulseHighFluxOffset+ModPulseLengthHighF/2.0) ;
    auto chopFrameOverlap1PhaseOffset=  (chopFrameOverlap1Offset)*14.0*360.0;

    auto chopFrameOverlap2Open= 1.0/14.0/InstLength*(chopFrameOverlap2Pos)*1.65 ;
    auto chopFrameOverlap2Offset=(  ( (chopFrameOverlap2Pos)/v_1+(chopFrameOverlap2Pos)/v_0)/2.0+PulseHighFluxOffset+ModPulseLengthHighF/2.0) ;
    auto chopFrameOverlap2PhaseOffset=  (chopFrameOverlap2Offset)*14.0*360.0;

/********** Bandwidth chopper calculations ****************/

//chopBW_t0= chopPulseOffset-chopPulseOpening/2.0 + (t_samp_0-(chopPulseOffset-chopPulseOpening/2.0)) / (InstLength-chopPulseDist) * (InstLength-chopBWPos) ;
//chopBW_t1= chopPulseOffset+chopPulseOpening/2.0 + (t_samp_1-(chopPulseOffset+chopPulseOpening/2.0)) / (InstLength-chopPulseDist) * (InstLength-chopBWPos);
    auto chopBW_t0= PulseHighFluxOffset+ModPulseLengthHighF/2.0 + chopBWPos/v_1;
    auto chopBW_t1=  PulseHighFluxOffset+ModPulseLengthHighF/2.0 + chopBWPos/v_0;

    auto chopBWOpen= 360.0/InstLength*(chopBWPos-chopPulseDist*1); //Here Jonas put a multiplier on the choppulsedist
    auto chopBWOffset=(chopBW_t0+chopBW_t1)/2.0;
    auto chopBWPhaseOffset=  (chopBWOffset)*14.0*360.0;

//    auto ps1 = std::make_pair(chopPulseFrequencyOrder * 14.0, chopPulsePhaseOffset);
//    auto ps2 = std::make_pair(chopPulseFrequencyOrder * 14.0, chopPulse2PhaseOffset);
//    auto fo1 = std::make_pair(14.0, chopFrameOverlap1PhaseOffset);
//    auto fo2 = std::make_pair(14.0, chopFrameOverlap2PhaseOffset);
//    // For some reason the instrument uses the delay time rather than the phase offset for the Bandwidth choppers
//    auto bw1 = std::make_pair(14.0, chopBWPhaseOffset);
//    auto bw2 = std::make_pair(-14.0, -chopBWPhaseOffset);

    std::map<std::string, chopper_parameters> cpm;
    cpm["ps1"] = {.speed=chopPulseFrequencyOrder * 14.0, .phase=chopPulsePhaseOffset, .angle=170.0, .path=chopPulseDist};//, .delay=chopPulseOffset};
    cpm["ps2"] = {.speed=chopPulseFrequencyOrder * 14.0, .phase=chopPulse2PhaseOffset, .angle=170.0, .path=chopPulseDist + 0.02}; //, .delay=chopPulseOffset};
    cpm["fo1"] = {.speed=14.0, .phase=chopFrameOverlap1PhaseOffset, .angle=38.26, .path=chopFrameOverlap1Pos}; // , .delay=chopFrameOverlap1Offset};
    cpm["fo2"] = {.speed=14.0, .phase=chopFrameOverlap2PhaseOffset, .angle=52.01, .path=chopFrameOverlap2Pos}; // .delay=chopFrameOverlap2Offset};
    cpm["bw1"] = {.speed=14.0, .phase=chopBWPhaseOffset, .angle=161.0, .path=chopBWPos}; // .delay=chopBWOffset};
    cpm["bw2"] = {.speed=-14.0, .phase=chopBWPhaseOffset, .angle=161.0, .path=chopBWPos + 0.02}; // .delay=chopBWOffset};

    return cpm;

//    std::map<std::string, double> choppers;
//    choppers["ps1speed"] = chopPulseFrequencyOrder * 14.0;
//    choppers["ps2speed"] = chopPulseFrequencyOrder * 14.0;
//    choppers["fo1speed"] = 14.0;
//    choppers["fo2speed"] = 14.0;
//    choppers["bw1speed"] = 14.0;
//    choppers["bw2speed"] = -14.0;
//    choppers["ps1phase"] = chopPulsePhaseOffset;
//    choppers["ps2phase"] = chopPulse2PhaseOffset;
//    choppers["fo1phase"] = chopFrameOverlap1PhaseOffset;
//    choppers["fo2phase"] = chopFrameOverlap2PhaseOffset;
//    choppers["bw1phase"] = chopBWPhaseOffset;
//    choppers["bw2phase"] = chopBWPhaseOffset;
//    choppers["ps1delay"] = chopPulseOffset;
//    choppers["ps2delay"] = chopPulseOffset;
//    choppers["fo1delay"] = chopFrameOverlap1Offset;
//    choppers["fo2delay"] = chopFrameOverlap2Offset;
//    choppers["bw1delay"] = chopBW_t0;
//    choppers["bw2delay"] = chopBW_t1;
//    return choppers;
}
