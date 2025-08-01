#pragma once
#include <vector>
#include <array>
#include <random>



struct config {

    double meltTemp; //K
    double dx; //um
    double dt; //s
    double kSolid; //W/m*k
    double kLiquid = kSolid * 0.45; //W/m*k
    int success;
    double startTemp; //K
    double cellWidth;//um
    double cellHeight;//um
    double cellArea = cellHeight * cellWidth;
    double particleVolFraction;
    double liqSolIntWidth; //um
    double grainIntWidth; //um

    int totalSteps;
    double liqSolIntE;
    double grainIntE;
    std::vector<int> steps;
    double phaseCoefficient;
    double barrierHeightPhase;
    double barrierHeightGrain;
    double phasePreCo = 0.75*liqSolIntE/(barrierHeightPhase*liqSolIntWidth);
    double grainPreCo = 0.75*grainIntE/(barrierHeightGrain*grainIntWidth);
    double particleSlowingCoefficient = 0;
    double phaseGradCo = 0.75*liqSolIntE*liqSolIntWidth;
    double basePlateTemp;

};




config inputConfig();


struct node {
    double temp;
    double phase;
    double particleComp;
    std::vector<double> grainPhases;
    
    node* neighbors[4]; //x-1,x+1,y-1,y+1
};

struct eulerAngles {
    double theta1;
    double phi;
    double theta2;
};
struct gridField {


    std::vector<std::vector<node>> grid;
    void buildGrid(int heightSteps, int widthSteps);
    void init(config modelConfig);
    
    int numGrains;
    std::vector<eulerAngles> orientations; // Theta1 Phi Theta2 Euler Angles
    node top;

    node bottom;


    void addGrain(std::vector<int> nucleus, config modelConf);

};





