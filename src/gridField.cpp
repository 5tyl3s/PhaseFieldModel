#include "gridField.hpp"
#include <iostream>
#include <random>
#include <cmath>

// Global grid dimensions - set at runtime
int GRID_ROWS = 500;
int GRID_COLS = 500;
int TOTAL_NODES = GRID_ROWS * GRID_COLS;

// This line actually creates the global object in memory
gridField globalField;
std::array<std::array<double,1000>,1000> tempGrid;

void gridField::buildGrid() {
    // Helper lambda to set neighbor indices
    auto setNeighbors = [this](int i, int j, int nodeIdx) {
        neighborIndices& nb = nodes.neighbors[nodeIdx];
        
        // left
        if (mConfig.periodicLeft) {
            int left_j = (j > 0) ? j - 1 : gridCols - 1;
            nb.idx[0] = i * gridCols + left_j;
        } else {
            if (j > 0) {
                nb.idx[0] = i * gridCols + (j - 1);
            } else {
                nb.idx[0] = totalNodes; // boundary marker
            }
        }
        
        // right
        if (mConfig.periodicRight) {
            int right_j = (j < gridCols - 1) ? j + 1 : 0;
            nb.idx[1] = i * gridCols + right_j;
        } else {
            if (j < gridCols - 1) {
                nb.idx[1] = i * gridCols + (j + 1);
            } else {
                nb.idx[1] = totalNodes; // boundary marker
            }
        }
        
        // up (top)
        if (mConfig.periodicTop) {
            int up_i = (i > 0) ? i - 1 : gridRows - 1;
            nb.idx[2] = up_i * gridCols + j;
        } else {
            if (i > 0) {
                nb.idx[2] = (i - 1) * gridCols + j;
            } else {
                nb.idx[2] = totalNodes; // boundary marker
            }
        }
        
        // down (bottom)
        if (mConfig.periodicBottom) {
            int down_i = (i < gridRows - 1) ? i + 1 : 0;
            nb.idx[3] = down_i * gridCols + j;
        } else {
            if (i < gridRows - 1) {
                nb.idx[3] = (i + 1) * gridCols + j;
            } else {
                nb.idx[3] = totalNodes + 1; // boundary marker
            }
        }
        
        // up-left
        int ul_i = (i > 0) ? i - 1 : (mConfig.periodicTop ? gridRows - 1 : gridRows - 1);
        int ul_j = (j > 0) ? j - 1 : (mConfig.periodicLeft ? gridCols - 1 : gridCols - 1);
        if ((i == 0 && !mConfig.periodicTop) || (j == 0 && !mConfig.periodicLeft)) {
            nb.idx[4] = totalNodes;
        } else {
            nb.idx[4] = ul_i * gridCols + ul_j;
        }
        
        // up-right
        int ur_i = (i > 0) ? i - 1 : (mConfig.periodicTop ? gridRows - 1 : gridRows - 1);
        int ur_j = (j < gridCols - 1) ? j + 1 : (mConfig.periodicRight ? 0 : 0);
        if ((i == 0 && !mConfig.periodicTop) || (j == gridCols - 1 && !mConfig.periodicRight)) {
            nb.idx[5] = totalNodes;
        } else {
            nb.idx[5] = ur_i * gridCols + ur_j;
        }
        
        // down-left
        int dl_i = (i < gridRows - 1) ? i + 1 : (mConfig.periodicBottom ? 0 : gridRows - 1);
        int dl_j = (j > 0) ? j - 1 : (mConfig.periodicLeft ? gridCols - 1 : gridCols - 1);
        if ((i == gridRows - 1 && !mConfig.periodicBottom) || (j == 0 && !mConfig.periodicLeft)) {
            nb.idx[6] = totalNodes + 1;
        } else {
            nb.idx[6] = dl_i * gridCols + dl_j;
        }
        
        // down-right
        int dr_i = (i < gridRows - 1) ? i + 1 : (mConfig.periodicBottom ? 0 : gridRows - 1);
        int dr_j = (j < gridCols - 1) ? j + 1 : (mConfig.periodicRight ? 0 : 0);
        if ((i == gridRows - 1 && !mConfig.periodicBottom) || (j == gridCols - 1 && !mConfig.periodicRight)) {
            nb.idx[7] = totalNodes + 1;
        } else {
            nb.idx[7] = dr_i * gridCols + dr_j;
        }
    };

    // Build all neighbor relationships
    for (int i = 0; i < gridRows; i++) {
        for (int j = 0; j < gridCols; j++) {
            int nodeIdx = i * gridCols + j;
            setNeighbors(i, j, nodeIdx);
        }
    }
}

void gridField::init(config modelConfig, int rows, int cols) {
    // Set runtime grid dimensions
    gridRows = rows;
    gridCols = cols;
    totalNodes = rows * cols;
    
    // Update global variables (for backward compatibility with code that references globals)
    GRID_ROWS = rows;
    GRID_COLS = cols;
    TOTAL_NODES = rows * cols;
    
    // Dynamically allocate all vectors to the required size
    nodes.resize(totalNodes);
    
    // store config for use in update()
    mConfig = modelConfig;
    canNucleate = true;  // Enable nucleation at start

    buildGrid();
    
    // Initialize boundary conditions
    top.grainPhases = {};
    top.phase = 0.00;
    top.exists = 0;
    top.activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    top.id = -1;
    
    bottom.grainPhases = {};
    bottom.phase = 1.00;
    bottom.exists = 0;
    bottom.particleComp = 0.0;
    bottom.activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    bottom.id = -1;
    
    numGrains = 0;
    
    static thread_local std::mt19937 rng{std::random_device{}()};
    static thread_local std::uniform_real_distribution<double> dist(0.5, 1.5);

    // Initialize all nodes
    for (int idx = 0; idx < totalNodes; idx++) {
        int i = idx / gridCols;
        int j = idx % gridCols;
        
        nodes.phase[idx] = 0.0;
        nodes.exists[idx] = 1;
        nodes.isDeactivated[idx] = 1;
        nodes.isParticleActive[idx] = 1;
        nodes.grainsHere[idx] = 0;
        nodes.grainPhases[idx] = {0,0,0,0,0,0,0,0,0};
        nodes.activeGrains[idx] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
        nodes.id[idx] = idx;
        nodes.heightPos[idx] = i;
        nodes.distFromCenter[idx] = 0.0f; // not used
        nodes.maxDist[idx] = 0.0f; // not used
        nodes.grainsToAdd[idx] = 0;
        nodes.hetNucleateHere[idx] = false;
        nodes.homoNucleateHere[idx] = false;
        nodes.sumGrains[idx] = 0;
        nodes.xPos[idx] = j;
        nodes.yPos[idx] = i;
        nodes.nucleatedHere[idx] = false;
        nodes.hasHetNucleated[idx] = false;
        nodes.hetNucleationSite[idx] = false;
        nodes.homoNucleationSite[idx] = false;
        nodes.maxGrainPhase[idx] = 0.0;
        nodes.particleComp[idx] = 0.0;
        if (modelConfig.hemisphericalCooling) {
            float centerCol = (gridCols - 1) / 2.0f;
            float maxRadius = std::sqrt((gridRows - 1)*(gridRows - 1) + centerCol*centerCol);
            float dist = std::sqrt(((gridRows - 1 - i) * (gridRows - 1 - i)) + ((j - centerCol) * (j - centerCol)));
            nodes.tempDistance[idx] = (dist < maxRadius) ? ((maxRadius - dist) * modelConfig.dx) : 0.0f;
        } else if (modelConfig.radialCooling) {
            nodes.tempDistance[idx] = (std::sqrt(std::pow(gridRows/2,2)+std::pow(gridCols/2,2))-1*std::sqrt(std::pow(i - gridRows/2, 2) + std::pow(j - gridCols/2, 2)) )* modelConfig.dx; // distance from center for radial cooling
        } else {
            nodes.tempDistance[idx] = ((gridRows - 1 - i) * modelConfig.dx); // vertical distance from top for planar cooling
        }
        nodes.baseTemp[idx] = static_cast<float>(modelConfig.startTemp + (nodes.tempDistance[idx] * modelConfig.tGrad));
        nodes.temp[idx] = nodes.baseTemp[idx];
    }
    
    // Initialize particles
    double particleVolume = (4.0/3.0) * 3.14159 * std::pow(modelConfig.particleRadius, 3);
    int numParts = static_cast<int>(std::round(gridRows * gridCols * modelConfig.particleVolFraction*modelConfig.dx*modelConfig.dx*modelConfig.thickness2D/particleVolume));
    std::cout << "Initializing with " << numParts << " particles.\n";
    for (int particle = 0; particle < numParts; particle++) {
        int x = rng() % gridCols;
        int y = rng() % gridRows;
        int idx = y * gridCols + x;
        nodes.particleComp[idx] += particleVolume/(modelConfig.dx*modelConfig.dx*modelConfig.thickness2D);
    }
}

void gridField::addGrain(int nodeIdx) {
    numGrains = numGrains + 1;
    std::cout << "There are now " << numGrains << " Grains\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(45,45);
    eulerAngles tempRots = {dist(gen),dist(gen),dist(gen)};
    
    // Initialize the grain at the nucleus location
    if (nodes.grainsHere[nodeIdx] < 9) {
        nodes.activeGrains[nodeIdx][nodes.grainsHere[nodeIdx]] = numGrains-1;
        nodes.grainPhases[nodeIdx][nodes.grainsHere[nodeIdx]] = 1;
        nodes.orientations[nodeIdx][nodes.grainsHere[nodeIdx]] = tempRots;
        nodes.grainsHere[nodeIdx]++;
    }
    nodes.phase[nodeIdx] = 1;

    // Add to neighbors
    for (int nb = 0; nb < 4; nb++) {
        int nbrIdx = nodes.neighbors[nodeIdx].idx[nb];
        
        // Skip boundaries and invalid indices
        if (nbrIdx >= TOTAL_NODES || nodes.exists[nbrIdx] == 0) continue;
        
        if (nodes.grainsHere[nbrIdx] < 9) {
            nodes.activeGrains[nbrIdx][nodes.grainsHere[nbrIdx]] = numGrains-1;
            nodes.orientations[nbrIdx][nodes.grainsHere[nbrIdx]] = tempRots;
            nodes.grainPhases[nbrIdx][nodes.grainsHere[nbrIdx]] = 1;
            nodes.grainsHere[nbrIdx]++;
        }
    }
}

void gridField::update(
    std::vector<float> &phaseDiffEn,
    std::vector<std::array<float,9>> &grainDiffEn,
    std::vector<float> &tempPartComp,
    std::vector<float> &tGrad,
    bool enableProfiling
) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    bool grainExists = 0;

    static long long totalFirstUpdateTime = 0;
    static long long totalParticleUpdateTime = 0;
    static long long totalPropagateGrainTime = 0;
    static int profileCount = 0;

    // Update temperature for all nodes
    #pragma omp parallel for
    for (int idx = 0; idx < totalNodes; idx++) {
        nodes.temp[idx] = tGrad[idx];
    }

    // Activate nodes when temperature drops below melt temp
    #pragma omp parallel for
    for (int idx = 0; idx < totalNodes; idx++) {
        if (nodes.temp[idx] < mConfig.meltTemp && nodes.isDeactivated[idx]) {
            nodes.isDeactivated[idx] = 0;
        }
    }

    // First pass: update phase and grain phases
    #pragma omp parallel for
    for (int idx = 0; idx < totalNodes; idx++) {
        if (nodes.isDeactivated[idx]) continue;
        nodes.phase[idx] = nodes.phase[idx] - pow((1-nodes.particleComp[idx]),3)*(mConfig.dt / pow(mConfig.dx,2)) * phaseDiffEn[idx]*4e-14;
        double grainHereCount = 0.0;
        float maxGrainPhase = 0.0f;

        for (int g = 0; g < nodes.grainsHere[idx]; g++) {
            nodes.grainPhases[idx][g] = nodes.grainPhases[idx][g] - (mConfig.dt / pow(mConfig.dx,2) * grainDiffEn[idx][g]*2e-14);
            if (nodes.grainPhases[idx][g] < 0.0) nodes.grainPhases[idx][g] = 0.0;
            if (nodes.grainPhases[idx][g] > 1.0) nodes.grainPhases[idx][g] = 1.0;
            grainHereCount = grainHereCount + pow(nodes.grainPhases[idx][g],2);
            maxGrainPhase = std::max(maxGrainPhase, nodes.grainPhases[idx][g]);
        }

        // Remove non-existing grain entries: any grain slot with zero phase or invalid id
        int writeIdx = 0;
        for (int readIdx = 0; readIdx < nodes.grainsHere[idx]; readIdx++) {
            bool isAlive = (nodes.grainPhases[idx][readIdx] > 0.0f) && (nodes.activeGrains[idx][readIdx] >= 0);
            if (isAlive) {
                if (writeIdx != readIdx) {
                    nodes.activeGrains[idx][writeIdx] = nodes.activeGrains[idx][readIdx];
                    nodes.grainPhases[idx][writeIdx] = nodes.grainPhases[idx][readIdx];
                    nodes.orientations[idx][writeIdx] = nodes.orientations[idx][readIdx];
                }
                writeIdx++;
            }
        }
        for (int clearIdx = writeIdx; clearIdx < 9; clearIdx++) {
            nodes.activeGrains[idx][clearIdx] = -1;
            nodes.grainPhases[idx][clearIdx] = 0.0f;
        }
        nodes.grainsHere[idx] = writeIdx;
        nodes.sumGrains[idx] = grainHereCount;
        nodes.maxGrainPhase[idx] = maxGrainPhase;
    }

    // Deactivate particle updates for solidified nodes
    #pragma omp parallel for
    for (int idx = 0; idx < totalNodes; idx++) {
        if (nodes.phase[idx] > 0.95) {
            nodes.isParticleActive[idx] = 0;
        }
    }

    // Deactivate nodes that are fully solidified with grains
    #pragma omp parallel for
    for (int idx = 0; idx < totalNodes; idx++) {
        if (nodes.phase[idx] > 0.95 && nodes.grainsHere[idx] > 0) {
            nodes.isDeactivated[idx] = 1;
        }
    }

    auto t1_start = std::chrono::steady_clock::now();

    const double particleMobilityLiquid = 2e-3;
    const double particleMobilitySolid = 1e-18;
    double dx = mConfig.dx;
    double denom = (dx * dx);
    
    std::vector<double> particleCompChange(totalNodes, 0.0);
    
    // Particle composition update
    #pragma omp parallel for
    for (int ptr = 0; ptr < totalNodes; ptr++) {
        if (!nodes.isParticleActive[ptr]) continue;
        double mu = tempPartComp[ptr];
        double sumNeighMu = 0.0;
        int nExist = 0;
        double selfM = (nodes.phase[ptr] * particleMobilitySolid) + ((1.0 - nodes.phase[ptr]) * particleMobilityLiquid);
        double totFlow = 0.0;
        double muFlux = 0.0;
        
        for (int nb = 0; nb < 4; nb++) {
            int nbrIdx = nodes.neighbors[ptr].idx[nb];
            double mu_nb = mu;
            double mNeigh = selfM;

            if (nbrIdx < totalNodes && nodes.exists[nbrIdx] != 0) {
                mu_nb = tempPartComp[nbrIdx];
                mNeigh = 0.5 * (
                    (nodes.phase[nbrIdx] * particleMobilitySolid + (1.0 - nodes.phase[nbrIdx]) * particleMobilityLiquid)
                    + selfM
                );
            }

            muFlux = mNeigh * (mu_nb - mu) / (mConfig.thickness2D * pow(dx,2));
            totFlow += muFlux;
        }
        particleCompChange[ptr] = 1*mConfig.dt * totFlow;
    }
    
    auto t1_end = std::chrono::steady_clock::now();
    if (enableProfiling) {
        totalFirstUpdateTime += std::chrono::duration_cast<std::chrono::microseconds>(t1_end - t1_start).count();
    }

    auto t2_start = std::chrono::steady_clock::now();
    
    // Apply particle composition changes
    #pragma omp parallel for
    for (int ptr = 0; ptr < totalNodes; ptr++) {
        if (!nodes.isParticleActive[ptr]) continue;
        nodes.particleComp[ptr] = nodes.particleComp[ptr] + particleCompChange[ptr];
    }

    auto t2_end = std::chrono::steady_clock::now();
    if (enableProfiling) {
        totalParticleUpdateTime += std::chrono::duration_cast<std::chrono::microseconds>(t2_end - t2_start).count();
    }

    auto t3_start = std::chrono::steady_clock::now();

    // Second pass: propagate grains and nucleation
    #pragma omp parallel for
    for (int ptr = 0; ptr < totalNodes; ptr++) {
        if (nodes.isDeactivated[ptr]) continue;
        if (nodes.phase[ptr] > 1.0) nodes.phase[ptr] = 1.0;
        if (nodes.phase[ptr] < 0.0) nodes.phase[ptr] = 0.0;
        
        for (int nb = 0; nb < 8; nb++) {
            int nbrIdx = nodes.neighbors[ptr].idx[nb];
            if (nbrIdx >= totalNodes || nodes.exists[nbrIdx] == 0) continue;
            
            double maxPhaseFromNeighbor = 0.0;
            for (int rg = 0; rg < nodes.grainsHere[nbrIdx]; rg++) {
                if (nodes.grainPhases[nbrIdx][rg] > maxPhaseFromNeighbor) {
                    maxPhaseFromNeighbor = nodes.grainPhases[nbrIdx][rg];
                }
            }
            
            if (maxPhaseFromNeighbor > 0.4) {
                for (int rg = 0; rg < nodes.grainsHere[nbrIdx]; rg++) {
                    if (nodes.grainPhases[nbrIdx][rg] >= maxPhaseFromNeighbor) {
                        int grainToAdd = nodes.activeGrains[nbrIdx][rg];
                        
                        // Check if we already have this grain
                        bool alreadyHere = false;
                        for (int lg = 0; lg < nodes.grainsHere[ptr]; lg++) {
                            if (nodes.activeGrains[ptr][lg] == grainToAdd) {
                                alreadyHere = true;
                                break;
                            }
                        }
                        
                        if (!alreadyHere && nodes.grainsToAdd[ptr] < 9) {
                            nodes.addGrainsHere[ptr][nodes.grainsToAdd[ptr]] = grainToAdd;
                            nodes.addGrainsOrientations[ptr][nodes.grainsToAdd[ptr]] = nodes.orientations[nbrIdx][rg];
                            nodes.grainsToAdd[ptr]++;
                        }
                    }
                }
            }
        }
        
        if (nodes.temp[ptr] < (mConfig.meltTemp)) {
            grainExists = 0;
            if (nodes.grainsHere[ptr] > 0 ) grainExists = 1;
            if (nodes.grainsToAdd[ptr] > 0) grainExists = 1;
            if (!grainExists && prob_dist(gen) < (1 - exp(pow(mConfig.dx,2)*mConfig.thickness2D * mConfig.dt * calcNucRate(nodes.temp[ptr], mConfig) * -1)) && nodes.phase[ptr] < 0.1) { 
                nodes.homoNucleateHere[ptr] = true;
            }
            if (!grainExists &&  nodes.temp[ptr] < mConfig.meltTemp-(mConfig.hetNucUnderCooling)&& prob_dist(gen) < nodes.particleComp[ptr]*mConfig.dt && !nodes.hasHetNucleated[ptr] && nodes.phase[ptr] < 0.1) { 
                nodes.hetNucleateHere[ptr] = true;
                std::cout << "Node " << ptr << " eligible for heterogeneous nucleation. Particle Comp: " << nodes.particleComp[ptr] << std::endl;
            }
        }
        
    }
    
    // Add new grains from neighbors and handle nucleation
    for (int ptr = 0; ptr < totalNodes; ptr++) {
        if (nodes.isDeactivated[ptr]) continue;
        double globAvailPart = nodes.particleComp[ptr]*(1-nodes.phase[ptr]);
        
        // Add grains from neighbors
        for (int ag = 0; ag < nodes.grainsToAdd[ptr]; ag++) {
            if (nodes.grainsHere[ptr] < 9) {
                int newGrainID = nodes.addGrainsHere[ptr][ag];
                nodes.activeGrains[ptr][nodes.grainsHere[ptr]] = newGrainID;
                nodes.orientations[ptr][nodes.grainsHere[ptr]] = nodes.addGrainsOrientations[ptr][ag];
                nodes.grainsHere[ptr]++;
            }
        }
        nodes.grainsToAdd[ptr] = 0;
        nodes.addGrainsHere[ptr] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
        
        if (nodes.hetNucleateHere[ptr]) {
            nodes.hetNucleateHere[ptr] = false;
            nodes.nucleatedHere[ptr] = true;
            nodes.hasHetNucleated[ptr] = true;
            nodes.hetNucleationSite[ptr] = true;
            addGrain(ptr);
            std::cout << "Heterogeneous Nucleation at Node " << ptr << " (x=" << nodes.xPos[ptr] << ", y=" << nodes.yPos[ptr] << ")\n";
        }
        if (nodes.homoNucleateHere[ptr]) {  
            nodes.nucleatedHere[ptr] = true;
            nodes.homoNucleationSite[ptr] = true;
            nodes.homoNucleateHere[ptr] = false;
            addGrain(ptr);
            std::cout << "Homogeneous Nucleation at Node " << ptr << " (x=" << nodes.xPos[ptr] << ", y=" << nodes.yPos[ptr] << ")\n";
        }
    }
    
    auto t3_end = std::chrono::steady_clock::now();
    if (enableProfiling) {
        totalPropagateGrainTime += std::chrono::duration_cast<std::chrono::microseconds>(t3_end - t3_start).count();
    }

    // Report performance
    if (enableProfiling) {
        profileCount++;
        if (profileCount >= 1000) {
            long long avgDiffEn = diffEnergyTime;
            long long avgFirstUpd = totalFirstUpdateTime / profileCount;
            long long avgPartUpd = totalParticleUpdateTime / profileCount;
            long long avgPropGrain = totalPropagateGrainTime / profileCount;
            
            std::cout << "\n=== Performance Profile (average over " << profileCount << " updates) ===" << std::endl;
            std::cout << "  Diff Energy Calculation: " << avgDiffEn << " μs" << std::endl;
            std::cout << "  First Update (Phase/Grain): " << avgFirstUpd << " μs" << std::endl;
            std::cout << "  Particle Composition Update: " << avgPartUpd << " μs" << std::endl;
            std::cout << "  Propagate Active Grain/Nucleation: " << avgPropGrain << " μs" << std::endl;
            std::cout << "  Total: " << (avgDiffEn + avgFirstUpd + avgPartUpd + avgPropGrain) << " μs" << std::endl;
            std::cout << "======================================\n" << std::endl;
            
            totalFirstUpdateTime = 0;
            totalParticleUpdateTime = 0;
            totalPropagateGrainTime = 0;
            profileCount = 0;
        }
    }
}

void gridField::recordDiffEnergyTime(long long timeMs) {
    diffEnergyTime = timeMs;
}

// Calculate nucleation rate based on temperature
float calcNucRate(double temp, config modelConf) {
    double k  = 1.380649e-23;        // Boltzmann constant (J/K)
    double h  = 6.62607015e-34;      // Planck constant (J·s)
    double NA = 6.02214076e23;       // Avogadro (1/mol)

    double dForceMo = (modelConf.drivingForceSlopek * temp
                     + modelConf.drivingForceIntercept)*5000/modelConf.molarVolume;

    double cellVolume = pow(modelConf.dx,2)*modelConf.thickness2D;
    double atomicDensity = NA / modelConf.molarVolume;
    double numAtoms = atomicDensity * cellVolume;
    double gamma = 0.05;
    
    double dGstar = (16.0 * 3.141592 * pow(gamma, 3.0)) /
                    (3.0 * dForceMo * dForceMo);
    
    double diffTerm = exp(-modelConf.diffusionActivationEnergy /
                          (8.314462618 * temp));
    
    double prefactor = (k * temp / h) * numAtoms;
    double rate = prefactor * diffTerm * exp(-dGstar / (k * temp));
    return (float)rate;
}
