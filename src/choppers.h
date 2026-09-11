//
// Created by gst on 23/10/23.
//

#ifndef BIFROST_CHOPPERS_CHOPPERS_H
#define BIFROST_CHOPPERS_CHOPPERS_H

#include <map>
#include <string>
#include "chopper.h"

auto bifrost(double E_0, double L_0, double chopPulseOpening) -> std::map<std::string, chopcal::Chopper>;

#endif //BIFROST_CHOPPERS_CHOPPERS_H
