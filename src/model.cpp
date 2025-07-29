#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>

void say_hello() {
    std::cout << "Hello From Model"<< std::endl;
}

void calcFreeEnergy(gridField field,config modelConfig) {

    for (int i = 1; i < modelConfig.steps[0]-1; i++) {
        for (int j = 1; j < modelConfig.steps[1]-1; j++) {
           // double fLocal = calcLocalFree(field.grid[i][j],modelConfig);
            //double fPhaseGrad = 0.5*modelConfig.phaseCoefficient*(field.grid[i][j].neighbors[0]->phase +field.grid[i][j].neighbors[1]->phase +field.grid[i][j].neighbors[2]->phase +field.grid[i][j].neighbors[3]->phase -4*field.grid[i][j].phase);//Von Nueman laplacian
            //double fGrainGrad = 0.5*static_cast<double>(calcGrainBoundaryEnergy(field.orientation,field.grid[i][j].grainPhases))*modelConfig.grainIntWidth;

        };
    };


};

void calcLocalFree(node nodeLoc, config mConfig) {
    //double energyPhase = mConfig.phasePreCo*(((1-nodeLoc.phase)*(1-nodeLoc.phase)*(ifLiq(nodeLoc.temp,mConfig.meltTemp)))+(nodeLoc.phase*nodeLoc.phase*(1-ifLiq(nodeLoc.temp,mConfig.meltTemp))));


};
void ifLiq(double temp, double meltTemp) {
    //0 Above melt, 1 below
   // double liqVal = 0.5*(1-(tanh((temp/meltTemp)-1)));
    //return liqVal;

};

std::array<std::array<double,3>,2> calcGrainBoundaryEnergy(eulerAngles orient, std::vector<double> grainPhase, std::array<double,2> gradient) {
    //Get grain <001> vector from euler angles
    std::array<double,3> directiion001 = {0,0,1};//BD
    std::array<double,3> direction010 = {0,1,0};//TD
    std::array<double, 3> newD001 = {sin(orient.theta1)*sin(orient.phi),-1*cos(orient.theta1)*sin(orient.phi),cos(orient.phi)};
    std::array<double,3> newD010 = {(-cos(orient.theta1)*sin(orient.theta2))-(sin(orient.theta1)*cos(orient.theta2)*cos(orient.phi)),((-1*sin(orient.theta1)*sin(orient.theta2))+(cos(orient.theta1)*cos(orient.theta2)*cos(orient.phi))),(cos(orient.theta2)*sin(orient.phi))};
    std::array<double, 2> gradNormal = {-1*gradient[1],1*gradient[0]};


    double en110 = 2.783;//J/m*m
    double en111 = 2.962;
    double en100 = 3.192;
    std::array<std::array<double,3>,2> out = {newD001,newD010};
    return out;



};

