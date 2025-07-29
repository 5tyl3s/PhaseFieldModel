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
    std::array<int,2> steps;
    double phaseCoefficient;
    double barrierHeightPhase;
    double barrierHeightGrain;
    double phasePreCo = 0.75*liqSolIntE/(barrierHeightPhase*liqSolIntWidth);
    double grainPreCo = 0.75*grainIntE/(barrierHeightGrain*grainIntWidth);
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
    
    std::vector<int> grains;
    eulerAngles orientation; // Theta1 Phi Theta2 Euler Angles
    std::random_device rd;

    static gridField addGrain(std::array<int,2> nucleus);
       
    
 
    
};





