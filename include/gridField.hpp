#pragma once
#include <array>
#include <vector>
#include "importConfig.hpp"

// Grid dimensions - declared as extern so they can be set at runtime
extern int GRID_ROWS;
extern int GRID_COLS;
extern int TOTAL_NODES;

struct eulerAngles {
    double theta1;
    double phi;
    double theta2;
};

// Helper struct for neighbor indices
struct neighborIndices {
    int idx[8];  // left, right, up, down, up-left, up-right, down-left, down-right
};

// Structure of Arrays for node data - improved cache locality and vectorization
// Uses dynamic allocation for scalability to multi-million node grids
struct nodeArrays {
    // Scalar fields - using float for temp/phase to save memory (50% reduction)
    std::vector<int> id;
    std::vector<float> temp;           // float instead of double (OK for thermal simulations)
    std::vector<float> phase;          // float is sufficient for phase field [0,1]
    std::vector<float> particleComp;
    std::vector<uint8_t> exists;       // uint8_t instead of double (8x smaller)
    std::vector<float> maxGrainPhase;
    std::vector<uint8_t> sumGrains;    // uint8_t sufficient for sum of grain phases
    std::vector<uint8_t> grainsHere;   // uint8_t sufficient (max 9 grains per node)
    std::vector<int> heightPos;
    std::vector<int> xPos;
    std::vector<int> yPos;
    
    // Boolean flags - using uint8_t instead of bool (more cache efficient)
    std::vector<uint8_t> hetNucleateHere;
    std::vector<uint8_t> homoNucleateHere;
    std::vector<uint8_t> nucleatedHere;
    std::vector<uint8_t> hasHetNucleated;
    std::vector<uint8_t> hetNucleationSite;
    std::vector<uint8_t> homoNucleationSite;
    
    // Per-grain arrays (9 grains per node)
    std::vector<std::array<int, 9>> addGrainsHere;         // int for grain IDs (can be -1 sentinel)
    std::vector<std::array<eulerAngles, 9>> addGrainsOrientations;
    std::vector<uint8_t> grainsToAdd;                          // uint8_t sufficient
    std::vector<std::array<eulerAngles, 9>> orientations;
    std::vector<std::array<float, 9>> grainPhases;             // float instead of double
    std::vector<std::array<int, 9>> activeGrains;
    
    // Neighbor indices (Moore neighborhood: 8 neighbors)
    std::vector<neighborIndices> neighbors;
    
    // Resize all vectors to match grid size
    void resize(int totalNodes) {
        id.resize(totalNodes);
        temp.resize(totalNodes);
        phase.resize(totalNodes);
        particleComp.resize(totalNodes);
        exists.resize(totalNodes);
        maxGrainPhase.resize(totalNodes);
        sumGrains.resize(totalNodes);
        grainsHere.resize(totalNodes);
        heightPos.resize(totalNodes);
        xPos.resize(totalNodes);
        yPos.resize(totalNodes);
        hetNucleateHere.resize(totalNodes);
        homoNucleateHere.resize(totalNodes);
        nucleatedHere.resize(totalNodes);
        hasHetNucleated.resize(totalNodes);
        hetNucleationSite.resize(totalNodes);
        homoNucleationSite.resize(totalNodes);
        addGrainsHere.resize(totalNodes);
        addGrainsOrientations.resize(totalNodes);
        grainsToAdd.resize(totalNodes);
        orientations.resize(totalNodes);
        grainPhases.resize(totalNodes);
        activeGrains.resize(totalNodes);
        neighbors.resize(totalNodes);
    }
};

// Helper function to calculate nucleation rate
float calcNucRate(double temp, config modelConf);

struct gridField {
    config mConfig;
    nodeArrays nodes;
    int numGrains;
    
    // Grid dimensions
    int gridRows, gridCols, totalNodes;
    
    // Boundary node data
    struct {
        std::array<int, 9> activeGrains;
        std::array<float, 9> grainPhases;  // float instead of double
        float phase;                        // float instead of double
        float particleComp;                 // float instead of double
        uint8_t exists;
        int id;
    } top, bottom;
    
    std::array<int,2> grainLocTemporary;
    long long diffEnergyTime;  // Store diff energy calculation time
    bool canNucleate;  // Flag to enable/disable nucleation when particles available

    void buildGrid();
    void init(config modelConfig, int rows = 500, int cols = 500);
    void addGrain(int nodeIdx);
    void update(
        std::vector<float> &phaseDiffEn,
        std::vector<std::array<float, 9>> &grainDiffEn,
        std::vector<float> &tempPartComp,
        std::vector<float> &tGrad,
        bool enableProfiling = true
    );
    void recordDiffEnergyTime(long long timeMs);
};

// Declare global variable (but don't allocate it here)
extern gridField globalField;