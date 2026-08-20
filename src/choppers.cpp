//
// Created by gst on 23/10/23.
//
#include <cmath>
#include <cstdio>
#include <vector>
#include <chopper-lib.h>
#include "choppers.h"

// The chopper structure is a run of doubles, so a library that means something else by
// its second field would accept these settings and quietly answer a different question.
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 20000
#error "chopcal sets chopper delays in seconds; chopper-lib 2.0.0 or newer is required"
#endif

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
    auto chopPulseSpeed = chopPulseFrequencyOrder * 14.0;
    auto chopPulseOffset=(chopPulseDist/v_1+chopPulseDist/v_0)/2.0+ModPulseLengthHighF/2.0+PulseHighFluxOffset;

    // The two disks turn together -- same speed, same direction -- and each is open for
    // as long as its slit takes to cross the beam. Running one behind the other leaves
    // only the part of that crossing they are both open for, so the burst is as short as
    // the lag makes it: a slit crossing less the lag. Lag them by the crossing less the
    // requested opening and the burst is the requested opening, at any speed.
    //
    // A disk is set by when its own slit centre is on the beam, so the two settings sit
    // half a burst either side of the band's arrival, each moved out by half a crossing.
    double chopPulseDelay = 0.0, chopPulse2Delay = 0.0;
    if (chopPulseFrequencyOrder == 0) {
        printf(" \n \n Warning: Pulse shaping chopper parked! Setting the offsets to zero");
    } else {
        auto chopPulseOpenTime = 170.0 / 360.0 / chopPulseSpeed;  // one slit, in seconds
        chopPulseDelay = chopPulseOffset + chopPulseOpening / 2.0 - chopPulseOpenTime / 2.0;
        chopPulse2Delay = chopPulseDelay - chopPulseOpening + chopPulseOpenTime;
    }


/*********** Frame Overlap chopper calculations ******************/
    // The offsets are already the times a band's middle reaches each disk, which is what
    // a chopper delay is. They were multiplied up by 14 * 360 to state them as an angle.
    auto chopFrameOverlap1Open= 1.0/14.0/InstLength*(chopFrameOverlap1Pos)*1.5 ;
    auto chopFrameOverlap1Offset=(  ( (chopFrameOverlap1Pos)/v_1+(chopFrameOverlap1Pos)/v_0)/2.0+PulseHighFluxOffset+ModPulseLengthHighF/2.0) ;

    auto chopFrameOverlap2Open= 1.0/14.0/InstLength*(chopFrameOverlap2Pos)*1.65 ;
    auto chopFrameOverlap2Offset=(  ( (chopFrameOverlap2Pos)/v_1+(chopFrameOverlap2Pos)/v_0)/2.0+PulseHighFluxOffset+ModPulseLengthHighF/2.0) ;

/********** Bandwidth chopper calculations ****************/

//chopBW_t0= chopPulseOffset-chopPulseOpening/2.0 + (t_samp_0-(chopPulseOffset-chopPulseOpening/2.0)) / (InstLength-chopPulseDist) * (InstLength-chopBWPos) ;
//chopBW_t1= chopPulseOffset+chopPulseOpening/2.0 + (t_samp_1-(chopPulseOffset+chopPulseOpening/2.0)) / (InstLength-chopPulseDist) * (InstLength-chopBWPos);
    auto chopBW_t0= PulseHighFluxOffset+ModPulseLengthHighF/2.0 + chopBWPos/v_1;
    auto chopBW_t1=  PulseHighFluxOffset+ModPulseLengthHighF/2.0 + chopBWPos/v_0;

    auto chopBWOpen= 360.0/InstLength*(chopBWPos-chopPulseDist*1); //Here Jonas put a multiplier on the choppulsedist
    // the middle of the band the bandwidth choppers pass, so their delay outright
    auto chopBWOffset=(chopBW_t0+chopBW_t1)/2.0;

    // The second bandwidth disk counter-rotates, and takes the same delay as the first:
    // both are on the beam at the middle of the band, whichever way they turn to get
    // there. Stated as a phase this needed the reader to know that chopper-lib divided
    // by the magnitude of the speed, so that the same positive angle meant the same
    // positive time for either sign.
    std::map<std::string, chopper_parameters> cpm;
    cpm["ps1"] = {.speed=chopPulseSpeed, .delay=chopPulseDelay, .angle=170.0, .path=chopPulseDist};
    cpm["ps2"] = {.speed=chopPulseSpeed, .delay=chopPulse2Delay, .angle=170.0, .path=chopPulseDist + 0.02};
    cpm["fo1"] = {.speed=14.0, .delay=chopFrameOverlap1Offset, .angle=38.26, .path=chopFrameOverlap1Pos};
    cpm["fo2"] = {.speed=14.0, .delay=chopFrameOverlap2Offset, .angle=52.01, .path=chopFrameOverlap2Pos};
    cpm["bw1"] = {.speed=14.0, .delay=chopBWOffset, .angle=161.0, .path=chopBWPos};
    cpm["bw2"] = {.speed=-14.0, .delay=chopBWOffset, .angle=161.0, .path=chopBWPos + 0.02};

    return cpm;
}
