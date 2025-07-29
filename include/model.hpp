#pragma once
#include "modelStartup.hpp"

std::array<std::array<double,3>,2> calcGrainBoundaryEnergy(eulerAngles orient, std::vector<double> grainPhase);