//
// Created by gst on 23/10/23.
//
#include <cmath>
#include <cstdio>
#include <vector>
#include <chopper-lib.h>
#include "choppers.h"
#include "constants.h"

// The chopper structure is a run of doubles, so a library that means something else by
// its second field would accept these settings and quietly answer a different question.
#if !defined(CHOPPER_LIB_VERSION) || CHOPPER_LIB_VERSION < 20000
#error "chopcal sets chopper delays in seconds; chopper-lib 2.0.0 or newer is required"
#endif

using namespace chopcal::constants;

auto bifrost(double E_0, double L_0, double chopPulseOpening)
-> std::map<std::string, chopper_parameters>
//-> std::map<std::string, std::map<std::string, double>>
//std::map<std::string, double>
{
// Transferred parameters
    double chopPulseFrequencyOrder=SOURCE_FREQUENCY; // Number of chopper pulses pr moderator pulse. It will automatically be reduced when nesesary and a warning will be written in the promt.
    double chopBWPos=BANDWIDTH_DISTANCE;  // Distance from pulse shapping choppers to BW Chopper

/*************************************** Chopper Variables  *******************************************/
    double chopFrameOverlap1Pos= FRAME_OVERLAP_1_DISTANCE;    // Distance from moderator to first frame owerlap chopper
    double chopFrameOverlap2Pos= FRAME_OVERLAP_2_DISTANCE;    // Distance from moderator to second frame owerlap chopper

// If there is set a value of L_0, overwrite E_0 and calculate E_0 from L_0
    if (L_0>0){
        E_0=H2_OVER_2M/(L_0*L_0);
    }


/************************************************/
/*                  Chopper calculations                    */
/************************************************/

    double PulseHighFluxOffset=PULSE_HIGH_FLUX_OFFSET; // Time from T0 to high pulse.
    double ModPulseLengthHighF=SOURCE_DURATION; // width of high pulse

    double InstLength=INSTRUMENT_LENGTH;
    double chopPulseDist= PULSE_SHAPING_DISTANCE;  // Distance fro moderator to Pulse chapping chopper

    /******* Check if pulse shapping chopper opening is large enough for requested frequency or reduce frequency *******/
    if  (chopPulseFrequencyOrder*chopPulseOpening > PULSE_SHAPING_ANGLE/DEGREES_PER_TURN/SOURCE_FREQUENCY) {
        chopPulseFrequencyOrder=floor(PULSE_SHAPING_ANGLE/DEGREES_PER_TURN/SOURCE_FREQUENCY/chopPulseOpening);
        printf(" \n \n Warning: Impossible combination of chopPulseFrequencyOrder and chopPulseOpening chosen, chopPulseFrequencyOrder reduced to: %f  \n", chopPulseFrequencyOrder);
    }

    auto lambda_1=sqrt(H2_OVER_2M/E_0);  /**** general chopper calculations **********/
    auto WavelengthBand = H_OVER_M/(InstLength-chopPulseDist)/SOURCE_FREQUENCY;
    auto lambda_0=lambda_1-WavelengthBand;
    auto v_0=H_OVER_M/lambda_1;
    auto v_1=H_OVER_M/lambda_0;

/***********  Pulse shaping chopper calculations **********/
    auto chopPulseSpeed = chopPulseFrequencyOrder * SOURCE_FREQUENCY;
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
        auto chopPulseOpenTime = PULSE_SHAPING_ANGLE / DEGREES_PER_TURN / chopPulseSpeed;  // one slit, in seconds
        chopPulseDelay = chopPulseOffset + chopPulseOpening / 2.0 - chopPulseOpenTime / 2.0;
        chopPulse2Delay = chopPulseDelay - chopPulseOpening + chopPulseOpenTime;
    }


/*********** Frame Overlap chopper calculations ******************/
    // The offsets are already the times a band's middle reaches each disk, which is what
    // a chopper delay is. They were multiplied up by 14 * 360 to state them as an angle.
    auto chopFrameOverlap1Offset=(  ( (chopFrameOverlap1Pos)/v_1+(chopFrameOverlap1Pos)/v_0)/2.0+PulseHighFluxOffset+ModPulseLengthHighF/2.0) ;

    auto chopFrameOverlap2Offset=(  ( (chopFrameOverlap2Pos)/v_1+(chopFrameOverlap2Pos)/v_0)/2.0+PulseHighFluxOffset+ModPulseLengthHighF/2.0) ;

/********** Bandwidth chopper calculations ****************/

//chopBW_t0= chopPulseOffset-chopPulseOpening/2.0 + (t_samp_0-(chopPulseOffset-chopPulseOpening/2.0)) / (InstLength-chopPulseDist) * (InstLength-chopBWPos) ;
//chopBW_t1= chopPulseOffset+chopPulseOpening/2.0 + (t_samp_1-(chopPulseOffset+chopPulseOpening/2.0)) / (InstLength-chopPulseDist) * (InstLength-chopBWPos);
    auto chopBW_t0= PulseHighFluxOffset+ModPulseLengthHighF/2.0 + chopBWPos/v_1;
    auto chopBW_t1=  PulseHighFluxOffset+ModPulseLengthHighF/2.0 + chopBWPos/v_0;

    // the middle of the band the bandwidth choppers pass, so their delay outright
    auto chopBWOffset=(chopBW_t0+chopBW_t1)/2.0;

    // The second bandwidth disk counter-rotates, and takes the same delay as the first:
    // both are on the beam at the middle of the band, whichever way they turn to get
    // there. Stated as a phase this needed the reader to know that chopper-lib divided
    // by the magnitude of the speed, so that the same positive angle meant the same
    // positive time for either sign.
    std::map<std::string, chopper_parameters> cpm;
    cpm["ps1"] = {.speed=chopPulseSpeed, .delay=chopPulseDelay, .angle=PULSE_SHAPING_ANGLE, .path=chopPulseDist};
    cpm["ps2"] = {.speed=chopPulseSpeed, .delay=chopPulse2Delay, .angle=PULSE_SHAPING_ANGLE, .path=chopPulseDist + PAIR_SEPARATION};
    cpm["fo1"] = {.speed=SOURCE_FREQUENCY, .delay=chopFrameOverlap1Offset, .angle=FRAME_OVERLAP_1_ANGLE, .path=chopFrameOverlap1Pos};
    cpm["fo2"] = {.speed=SOURCE_FREQUENCY, .delay=chopFrameOverlap2Offset, .angle=FRAME_OVERLAP_2_ANGLE, .path=chopFrameOverlap2Pos};
    cpm["bw1"] = {.speed=SOURCE_FREQUENCY, .delay=chopBWOffset, .angle=BANDWIDTH_ANGLE, .path=chopBWPos};
    cpm["bw2"] = {.speed=-SOURCE_FREQUENCY, .delay=chopBWOffset, .angle=BANDWIDTH_ANGLE, .path=chopBWPos + PAIR_SEPARATION};

    return cpm;
}
