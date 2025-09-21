#pragma once
#include <array>
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
    int heightPos;

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
    int numGrains;
    node top;
    node bottom;
    std::array<int,2> grainLocTemporary;

    void buildGrid();
    void init(config modelConfig);
    void addGrain(node* nucleus);
    void update(
        std::array<double,1000000> &phaseDiffEn,
        std::array<std::array<double,9>,1000000> &grainDiffEn,
        std::array<double,1000000> &tempPartComp,
        std::array<double,1000000> &tGrad
    );
};

// Declare global variable (but don't allocate it here)
extern gridField globalField;
