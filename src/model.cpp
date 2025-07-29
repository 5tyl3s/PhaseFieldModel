#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>

void say_hello() {
    std::cout << "Hello From Model"<< std::endl;
}

double calcFreeEnergy(gridField field,config modelConfig) {

    for (int i = 1; i < modelConfig.steps[0]-1; i++) {
        for (int j = 1; j < modelConfig.steps[1]-1; j++) {
            double fLocal = calcLocalFree(field.grid[i][j],modelConfig);
            double fPhaseGrad = 0.5*modelConfig.phaseCoefficient*(field.grid[i][j].neighbors[0]->phase +field.grid[i][j].neighbors[1]->phase +field.grid[i][j].neighbors[2]->phase +field.grid[i][j].neighbors[3]->phase -4*field.grid[i][j].phase);//Von Nueman laplacian
            double fGrainGrad = 0.5*calcGrainBoundaryEnergy(field.orientation,field.grid[i][j].grainPhases)*modelConfig.grainIntWidth;

        };
    };


};

double calcLocalFree(node nodeLoc, config mConfig) {
    double energyPhase = mConfig.phasePreCo*(((1-nodeLoc.phase)*(1-nodeLoc.phase)*(ifLiq(nodeLoc.temp,mConfig.meltTemp)))+(nodeLoc.phase*nodeLoc.phase*(1-ifLiq(nodeLoc.temp,mConfig.meltTemp))));


};
double ifLiq(double temp, double meltTemp) {
    //0 Above melt, 1 below
    double liqVal = 0.5*(1-(tanh((temp/meltTemp)-1)));
    return liqVal;

};

double calcGrainBoundaryEnergy(std::vector<std::vector<double>> orient, std::vector<double> grainPhase) {

};

