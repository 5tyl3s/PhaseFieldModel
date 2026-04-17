#pragma once
#include "modelStartup.hpp"
#include "gridField.hpp"
double underCool(double temp, double meltTemp);

double calcGrainBoundaryEnergy(eulerAngles orient, const std::array<double,3>& gradient);

std::vector<std::vector<double>> calcFreeEnergy(gridField& field, config modelConfig);

std::array<float,9> calcGrainDiffEnergy(int nodeIdx, config mConfig);

float calcPhaseDiffEnergy(int nodeIdx, config modelConfig);

double calcParticleCompDiff(int nodeIdx, config& modelConfig);

double calcTemp(int nodeIdx, config& modelConfig, int t);

void resetEnergyProfilingStats();
void printEnergyProfilingStats();

std::array<double,3> eulerRotate(eulerAngles orient, std::array<double,3> rotatedVector);

double dotAngle(std::array<double,3> vec1, std::array<double,3> vec2);
bool readCoeffs(const std::string &filename);


inline double getPotentialWithLaplacian(gridField& field, config& modelConfig, int i, int j);