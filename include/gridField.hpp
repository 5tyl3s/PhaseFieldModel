#pragma once
#include <array>
#include "importConfig.hpp"

// Define grid dimensions in one place. Change these values to the desired grid size.
constexpr int GRID_ROWS = 400;
constexpr int GRID_COLS = 400;
constexpr int TOTAL_NODES = GRID_ROWS * GRID_COLS;

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
    double maxGrainPhase;  // Track the maximum grain phase at this node for thresholding
    int sumGrains;
    int grainsHere;
    int heightPos;
    int xPos;
    int yPos;
    bool hetNucleateHere;
    bool homoNucleateHere;
    bool nucleatedHere;
    bool hasHetNucleated;  // Permanent flag: node has already undergone heterogeneous nucleation
    bool hetNucleationSite;  // Permanent tracking flag: this node was a het nucleation site (visualize red)
    bool homoNucleationSite;  // Permanent tracking flag: this node was a homo nucleation site (visualize pink)
    std::array<int,9> addGrainsHere;
    std::array<eulerAngles, 9> addGrainsOrientations;
    int grainsToAdd;

    std::array<eulerAngles, 9> orientations;
    std::array<double, 9> grainPhases;
    std::array<int, 9> activeGrains;

    // Use 8 neighbors (Moore): left, right, up, down, up-left, up-right, down-left, down-right
    std::array<node*, 8> neighbors;
    float calcNucRate(config modelConf);
    
};

struct gridField {
    config mConfig;
    std::array<std::array<node, GRID_COLS>, GRID_ROWS> grid;
    std::array<node*, TOTAL_NODES> allNodes;
    int numGrains;
    node top;
    node bottom;
    std::array<int,2> grainLocTemporary;
    long long diffEnergyTime;  // Store diff energy calculation time
    bool canNucleate;  // Flag to enable/disable nucleation when particles available

    void buildGrid();
    void init(config modelConfig);
    void addGrain(node* nucleus);
    void update(
        std::array<double, TOTAL_NODES> &phaseDiffEn,
        std::array<std::array<double,9>, TOTAL_NODES> &grainDiffEn,
        std::array<double, TOTAL_NODES> &tempPartComp,
        std::array<double, TOTAL_NODES> &tGrad,
        bool enableProfiling = true
    );
    void recordDiffEnergyTime(long long timeMs);};

// Declare global variable (but don't allocate it here)
extern gridField globalField;