#pragma once
#include "modelStartup.hpp"
double ifLiq(double temp, double meltTemp);
double calcGrainBoundaryEnergy(eulerAngles orient, std::vector<double> gradient);
double calcLocalFree(node nodeLoc, config mConfig);
std::vector<std::vector<double>> calcFreeEnergy(gridField field,config modelConfig);
double compOtherGrains(int notIndex,node Node, int numGrains);
std::vector<std::vector<std::vector<double>>> calcGrainDiffEnergy(gridField field, config modelConfig);
std::vector<std::vector<double>> calcPhaseDiffEnergy(gridField field, config modelConfig);
std::vector<double> eulerRotate(eulerAngles orient, std::vector<double> rotatedVector);
double dotAngle(std::vector<double> vec1, std::vector<double> vec2);
std::vector<std::vector<double>> calcTempDiff(gridField& field, const config& modelConfig);