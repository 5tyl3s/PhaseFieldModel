#pragma once
#include "modelStartup.hpp"
#include "gridField.hpp"
double ifLiq(double temp, double meltTemp);

double calcGrainBoundaryEnergy(eulerAngles orient, std::vector<double> gradient);

double calcLocalFree(node nodeLoc, config mConfig);

std::vector<std::vector<double>> calcFreeEnergy(gridField field,config modelConfig);

double compOtherGrains(int notIndex,node Node, int numGrains);

std::array<double,9> calcGrainDiffEnergy(node* nd, config mConfig);

double calcPhaseDiffEnergy(node* nd, config modelConfig);

std::array<double,3> eulerRotate(eulerAngles orient, std::array<double,3> rotatedVector);

double dotAngle(std::array<double,3> vec1, std::array<double,3> vec2);
double calcParticleCompDiff(node* nd, config& modelConfig);
double calcTemp(node* nd, config& modelConfig,int t);


inline double getPotentialWithLaplacian(gridField& field, config& modelConfig, int i, int j);