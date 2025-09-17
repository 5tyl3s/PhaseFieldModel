#pragma once
#include <vector>
#include <array>
#include <random>
#include "importConfig.hpp"







struct eulerAngles {
    double theta1;
    double phi;
    double theta2;
};




struct node {
    int id;
    double temp;
    double phase;
    double particleComp;
    double exists;
    int grainsHere;
    std::array<eulerAngles, 9> orientations;
    std::array<double, 9> grainPhases;
    std::array<int, 9> activeGrains;

    std::array<node*, 4> neighbors; //x-1,x+1,y-1,y+1
    float calcNucRate(config modelConf);
    
};



struct gridField {
    config mConfig;
    std::array<std::array<node,1000>,1000> grid;
    std::array<node*,1000000> allNodes;
    void buildGrid();
    void init(config modelConfig);
    

    int numGrains;
     // Theta1 Phi Theta2 Euler Angles
    node top;
    node bottom;
    std::array<int,2> grainLocTemporary;

    void addGrain(node* nucleus);

    

    // Updated update function for quadrant support:
    void update(
    std::array<double,1000000> &phaseDiffEn,
    std::array<double,1000000> &tempGrad,
    std::array<std::array<double,9>,1000000> &grainDiffEn,
    std::array<double,1000000> &tempPartComp
    );
};



std::array<double,1000000> updateTemp(double tGrad, double tRate, int timeStep, double dx, double dt, int iSteps, int jSteps, double startTemp, double minTemp);

