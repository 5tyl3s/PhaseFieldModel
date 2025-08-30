#pragma once
#include <vector>
#include <array>
#include <random>




struct config {

    double meltTemp; //K
    double dx; //um
    double dt; //s
    double kSolid; //W/m*k
    double kLiquid; //W/m*k
    int success;
    double startTemp; //K
    double cellWidth;//um
    double cellHeight;//um
    double cellArea;
    double particleVolFraction;
    double liqSolIntWidth; //um
    double grainIntWidth; //um
    int timeSteps;
    double underCoolReq;
    int totalSteps;
    double liqSolIntE;
    double grainIntE;
    std::vector<int> steps;
    double phaseCoefficient;
    double barrierHeightPhase;
    double barrierHeightGrain;
    double phasePreCo;
    double grainPreCo;
    double particleSlowingCoefficient;
    double phaseGradCo;
    double basePlateTemp;
    double heatCapacity;
    double density;
    double grainGradCo;
    double homoNucCoeff;
    double minTemp;
    double tGrad;
    double coolingRate;

};




config inputConfig();


struct node {
    int id;
    double temp;
    double phase;
    double particleComp;
    std::vector<double> grainPhases;
    std::vector<int> activeGrains;
    
    node* neighbors[4]; //x-1,x+1,y-1,y+1
};

struct eulerAngles {
    double theta1;
    double phi;
    double theta2;
};
struct gridField {
    std::vector<std::vector<node>> grid;
    std::vector<std::vector<node*>> grainActiveNodes;
    void buildGrid(int heightSteps, int widthSteps);
    void init(config modelConfig);
    

    int numGrains;
    std::vector<eulerAngles> orientations; // Theta1 Phi Theta2 Euler Angles
    node top;
    node bottom;
    std::vector<int> grainLocTemporary;

    void addGrain(std::vector<int> nucleus, config modelConf);

    // Updated update function for quadrant support:
    void update(
        std::vector<std::vector<double>> &phaseDiffEn,
        std::vector<std::vector<double>> &tempGrad,
        std::vector<std::vector<std::vector<double>>> &grainDiffEn,
        config &modelConf,
        int i_start, int i_end, int j_start, int j_end
    );
    void resizeGrainDiffEn(std::vector<std::vector<std::vector<double>>>& grainDiffEn, int numGrains) {
        for (auto& row : grainDiffEn) {
            for (auto& cell : row) {
                if (cell.size() != numGrains)
                    cell.resize(numGrains, 0.0);
            }
        }
    };
    
};




std::vector<std::vector<double>> updateTemp(double tGrad, double tRate, int timeStep, double dx, double dt, int iSteps, int jSteps, double startTemp, double minTemp);

