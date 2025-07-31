#pragma once
#include "modelStartup.hpp"
double ifLiq(double temp, double meltTemp);
double calcGrainBoundaryEnergy(eulerAngles orient, std::array<double,2> gradient);
double calcLocalFree(node nodeLoc, config mConfig);
std::vector<std::vector<double>> calcFreeEnergy(gridField field,config modelConfig);
double compOtherGrains(int notIndex,node Node, int numGrains);
std::vector<std::vector<std::vector<double>>> calcGrainDiffEnergy(gridField field, config modelConfig);
std::vector<std::vector<double>> calcPhaseDiffEnergy(gridField field, config modelConfig);
std::array<double,3> eulerRotate(eulerAngles orient, std::array<double,3> rotatedVector);
double dotAngle(std::array<double,3> vec1, std::array<double,3> vec2);
std::vector<std::vector<double>> calcTempDiff(gridField field, config modelConfig);