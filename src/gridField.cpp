#include "gridField.hpp"
#include <iostream>
#include <random>
#include <cmath>

// This line actually creates the global object in memory
gridField globalField;
std::array<std::array<double,1000>,1000> tempGrid;

void gridField::buildGrid() {

    // top row
    for (int j = 0; j < GRID_COLS; j++) {
        node& gp = grid[0][j];
        // left, right, up, down, up-left, up-right, down-left, down-right
        gp.neighbors[0] = (j > 0) ? &grid[0][j - 1] : &grid[0][GRID_COLS - 1]; // left
        gp.neighbors[1] = (j < GRID_COLS - 1) ? &grid[0][j + 1] : &grid[0][0]; // right
        gp.neighbors[2] = &top; // up
        gp.neighbors[3] = &grid[1][j]; // down
        gp.neighbors[4] = &top; // up-left (top row -> top)
        gp.neighbors[5] = &top; // up-right
        gp.neighbors[6] = &grid[1][ (j > 0) ? j - 1 : GRID_COLS - 1 ]; // down-left
        gp.neighbors[7] = &grid[1][ (j < GRID_COLS - 1) ? j + 1 : 0 ]; // down-right
         allNodes[j] = &gp;
         gp.heightPos = 0;
    }

    // interior rows
    for (int i = 1; i < GRID_ROWS - 1; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            node& gp = grid[i][j];
            // left, right, up, down, up-left, up-right, down-left, down-right
            gp.neighbors[0] = (j > 0) ? &grid[i][j - 1] : &grid[i][GRID_COLS - 1]; // left
            gp.neighbors[1] = (j < GRID_COLS - 1) ? &grid[i][j + 1] : &grid[i][0]; // right
            gp.neighbors[2] = &grid[i - 1][j]; // up
            gp.neighbors[3] = &grid[i + 1][j]; // down
            gp.neighbors[4] = &grid[i - 1][ (j > 0) ? j - 1 : GRID_COLS - 1 ]; // up-left
            gp.neighbors[5] = &grid[i - 1][ (j < GRID_COLS - 1) ? j + 1 : 0 ];       // up-right
            gp.neighbors[6] = &grid[i + 1][ (j > 0) ? j - 1 : GRID_COLS - 1 ]; // down-left
            gp.neighbors[7] = &grid[i + 1][ (j < GRID_COLS - 1) ? j + 1 : 0 ];       // down-right
             allNodes[i * GRID_COLS + j] = &gp;
             gp.heightPos = i;
        }
    }

    // bottom row
    for (int j = 0; j < GRID_COLS; j++) {
        node& gp = grid[GRID_ROWS - 1][j];
        gp.neighbors[0] = (j > 0) ? &grid[GRID_ROWS - 1][j - 1] : &grid[GRID_ROWS - 1][GRID_COLS - 1]; // left
        gp.neighbors[1] = (j < GRID_COLS - 1) ? &grid[GRID_ROWS - 1][j + 1] : &grid[GRID_ROWS - 1][0]; // right
        gp.neighbors[2] = &grid[GRID_ROWS - 2][j]; // up
        gp.neighbors[3] = &bottom; // down
        gp.neighbors[4] = &grid[GRID_ROWS - 2][ (j > 0) ? j - 1 : GRID_COLS - 1 ]; // up-left
        gp.neighbors[5] = &grid[GRID_ROWS - 2][ (j < GRID_COLS - 1) ? j + 1 : 0 ];       // up-right
        gp.neighbors[6] = &bottom; // down-left (bottom row -> bottom)
        gp.neighbors[7] = &bottom; // down-right
         allNodes[(GRID_ROWS - 1) * GRID_COLS + j] = &gp;
         gp.heightPos = (GRID_ROWS - 1);
    }

};

void gridField::init(config modelConfig) {

    // store config for use in update()
    mConfig = modelConfig;
    canNucleate = true;  // Enable nucleation at start

    buildGrid();
    top.grainPhases = {};
    top.phase = 0.00;
    top.exists = 0;
    top.particleComp = 0.0;  // No-flux boundary: initialize to safe value
    bottom.grainPhases = {};
    bottom.phase = 1.00;
    bottom.exists = 0;
    bottom.particleComp = 0.0;  // No-flux boundary: initialize to safe value
    // mark slots as empty with -1 (valid grain IDs start at 0)
    top.activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    bottom.activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    numGrains = 0;
    top.exists = 0; // no-flux boundary
    bottom.exists = 0; // no-flux boundary
    
    static thread_local std::mt19937 rng{std::random_device{}()};
    static thread_local std::uniform_real_distribution<double> dist(0.9, 1.1);



    for (int ik = 0; ik < GRID_ROWS; ik++) {
        for (int jk = 0; jk < GRID_COLS; jk++) {
            grid[ik][jk].phase = 0.0;
            grid[ik][jk].particleComp = modelConfig.particleVolFraction *dist(rng);
            grid[ik][jk].exists = 1;
            grid[ik][jk].grainsHere = 0;
            grid[ik][jk].grainPhases = {0,0,0,0,0,0,0,0,0};
            grid[ik][jk].activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
            grid[ik][jk].id = ik * GRID_COLS + jk;
            grid[ik][jk].heightPos = ik;
            grid[ik][jk].temp = 0.0;
            // Initialize temperature: startTemp + vertical temperature gradient
            grid[ik][jk].temp = modelConfig.startTemp + (ik * modelConfig.tGrad * modelConfig.dx);
            grid[ik][jk].grainsToAdd = 0;
            grid[ik][jk].hetNucleateHere = false;
            grid[ik][jk].homoNucleateHere = false;
            grid[ik][jk].sumGrains = 0;
            grid[ik][jk].xPos = jk;
            grid[ik][jk].yPos = ik;
            grid[ik][jk].nucleatedHere = false;
            grid[ik][jk].hasHetNucleated = false;
            grid[ik][jk].hetNucleationSite = false;
            grid[ik][jk].homoNucleationSite = false;

        }
    }
}

void gridField::addGrain(node* nucleus) {
    
    numGrains = numGrains + 1;
    std::cout << "There are now "<<numGrains << "Grains\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(45,45);
    eulerAngles tempRots = {dist(gen),dist(gen),dist(gen)};
    // Initialize the grain at the nucleus location
    if (nucleus->grainsHere < 9) {
        nucleus->activeGrains[nucleus->grainsHere] = numGrains-1;
        nucleus->grainPhases[nucleus->grainsHere] = 1;
        nucleus->orientations[nucleus->grainsHere] = tempRots;
        nucleus->grainsHere++;
    }
    nucleus->phase = 1;

    if (nucleus->neighbors[0]->grainsHere < 9) {
        nucleus->neighbors[0]->activeGrains[nucleus->neighbors[0]->grainsHere] = numGrains-1;
        nucleus->neighbors[0]->orientations[nucleus->neighbors[0]->grainsHere] = tempRots;
        nucleus->neighbors[0]->grainPhases[nucleus->neighbors[0]->grainsHere] = 1;
        nucleus->neighbors[0]->grainsHere++;
    }
    
    if (nucleus->neighbors[1]->grainsHere < 9) {
        nucleus->neighbors[1]->activeGrains[nucleus->neighbors[1]->grainsHere] = numGrains-1;
        nucleus->neighbors[1]->orientations[nucleus->neighbors[1]->grainsHere] = tempRots;
        nucleus->neighbors[1]->grainPhases[nucleus->neighbors[1]->grainsHere] = 1;
        nucleus->neighbors[1]->grainsHere++;
    }




    if (nucleus->neighbors[2]->exists != 0) {
        if (nucleus->neighbors[2]->grainsHere < 9) {
            nucleus->neighbors[2]->activeGrains[nucleus->neighbors[2]->grainsHere] = numGrains-1;
            nucleus->neighbors[2]->orientations[nucleus->neighbors[2]->grainsHere] = tempRots;
            nucleus->neighbors[2]->grainPhases[nucleus->neighbors[2]->grainsHere] = 1;
            nucleus->neighbors[2]->grainsHere++;
        }
    }
    if (nucleus->neighbors[3]->exists != 0) {
        if (nucleus->neighbors[3]->grainsHere < 9) {
            nucleus->neighbors[3]->activeGrains[nucleus->neighbors[3]->grainsHere] = numGrains-1;
            nucleus->neighbors[3]->orientations[nucleus->neighbors[3]->grainsHere] = tempRots;
            nucleus->neighbors[3]->grainPhases[nucleus->neighbors[3]->grainsHere] = 1;
            nucleus->neighbors[3]->grainsHere++;
        }
    }



}





void gridField::update(
    std::array<double, TOTAL_NODES> &phaseDiffEn,
    std::array<std::array<double,9>, TOTAL_NODES> &grainDiffEn,
    std::array<double, TOTAL_NODES> &tempPartComp, // now holds chemical potential μ at each node
    std::array<double, TOTAL_NODES> &tGrad,
    bool enableProfiling
) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    bool grainExists = 0;

    // Performance profiling
    static long long totalFirstUpdateTime = 0;
    static long long totalParticleUpdateTime = 0;
    static long long totalPropagateGrainTime = 0;
    static int profileCount = 0;
    

    // use compile-time TOTAL_NODES from header
    node* n;
    #pragma omp parallel for private(n)
    for (int node = 0; node < TOTAL_NODES; node++) {
        n = globalField.allNodes[node];
        // make sure pointer is valid (defensive)
        if (!n) continue;
        
        // update phase
        n->phase = n->phase -(mConfig.dt / pow(mConfig.dx,2)) * phaseDiffEn[node]*5e-15;
        double grainHereCount = 0.0;

        // update grain phases
        for (int g = 0; g < 9; g++) {
            n->grainPhases[g] = n->grainPhases[g] - (mConfig.dt / pow(mConfig.dx,2) * grainDiffEn[node][g]*5e-14);
            grainHereCount = grainHereCount + pow(n->grainPhases[g],2);
            

        }
        n->sumGrains = grainHereCount;
    }

    // Profiling: Track time for first update (phase and grain updates)
    auto t1_start = std::chrono::steady_clock::now();

    // THERMODYNAMIC particle update: dC/dt = mobility * laplacian(μ)
    // Particles flow DOWN the chemical potential gradient to minimize free energy
    // μ is computed in calcParticleCompDiff() and stored in tempPartComp array
    const double particleMobilityLiquid = 1e-6;   // high mobility in liquid
    const double particleMobilitySolid = 1e-18;    // lower mo6bility in solid
    double dx = mConfig.dx;
    double denom = (dx * dx);
    
    // Store composition changes
    std::array<double, TOTAL_NODES> particleCompChange;
    particleCompChange.fill(0.0);
    
    #pragma omp parallel for
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        if (!n) continue;

        // tempPartComp[ptr] contains the chemical potential μ at this node
        double mu = tempPartComp[ptr];

        // Compute 4-point Laplacian of chemical potential
        double sumNeighMu = 0.0;
        double sumNeighMobility = 0;
        int nExist = 0;
        double mNeigh;
        double selfM = (n->phase * particleMobilitySolid) + ((1.0 - n->phase) * particleMobilityLiquid);
        double totFlow = 0.0;
        for (int nb = 0; nb < 4; nb++) {  // ONLY 4-point: left, right, up, down
           node* nbr = n->neighbors[nb];
           if (!nbr || nbr->exists == 0) continue;
            mNeigh = 0.5*((nbr->phase * particleMobilitySolid) + ((1.0 - nbr->phase) * particleMobilityLiquid)+selfM);
            double muFlux = mNeigh*(tempPartComp[nbr->id] - mu)/pow(dx,2);
            totFlow += muFlux;
        }
        // Choose mobility based on phase (liquid vs solid)
       
        // Store change
        particleCompChange[ptr] = mConfig.dt * totFlow/(dx*dx);
    }
    
    auto t1_end = std::chrono::steady_clock::now();
    if (enableProfiling) {
        totalFirstUpdateTime += std::chrono::duration_cast<std::chrono::microseconds>(t1_end - t1_start).count();
    }

    // Profiling: Track time for particle composition updates
    auto t2_start = std::chrono::steady_clock::now();
    
    // Apply changes TOGETHER to preserve global sum
    double globAvailPart = 0.0;
    #pragma omp parallel for
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        if (!n) continue;
        double c = n->particleComp;
        n->particleComp = n->particleComp + (c*(1-c)) * particleCompChange[ptr];
        
        
        // Clamp to [0, 1] to prevent unphysical values

    }


    auto t2_end = std::chrono::steady_clock::now();
    if (enableProfiling) {
        totalParticleUpdateTime += std::chrono::duration_cast<std::chrono::microseconds>(t2_end - t2_start).count();
    }

    // Profiling: Track time for grain propagation and nucleation
    auto t3_start = std::chrono::steady_clock::now();

    // second pass: clamp values and propagate active grains / nucleation
    #pragma omp parallel for
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        if (!n) continue;
        if (n->phase > 1.0) n->phase = 1.0;
        if (n->phase < 0.0) n->phase = 0.0;
        for (int nb = 0; nb < 8; nb++) {
            node* nbr = n->neighbors[nb];
            if (!nbr || nbr->exists == 0) continue;
            
            // Find the highest grain phase from this neighbor to use as threshold
            double maxPhaseFromNeighbor = 0.0;
            for (int rg = 0; rg < nbr->grainsHere; rg++) {
                if (nbr->grainPhases[rg] > maxPhaseFromNeighbor) {
                    maxPhaseFromNeighbor = nbr->grainPhases[rg];
                }
            }
            
            // Add ALL grains with phase >= best phase threshold
            if (maxPhaseFromNeighbor > 0.4) {
                for (int rg = 0; rg < nbr->grainsHere; rg++) {
                    // Only propagate grains with phase >= best phase threshold
                    if (nbr->grainPhases[rg] >= maxPhaseFromNeighbor) {
                        int grainToAdd = nbr->activeGrains[rg];
                        
                        // Check if we already have this grain
                        bool alreadyHave = false;
                        for (int g = 0; g < n->grainsHere; g++) {
                            if (n->activeGrains[g] == grainToAdd) {
                                alreadyHave = true;
                                break;
                            }
                        }
                        
                        // Check if already in pending list
                        bool alreadyAdded = false;
                        for (int ag = 0; ag < n->grainsToAdd; ag++) {
                            if (n->addGrainsHere[ag] == grainToAdd) {
                                alreadyAdded = true;
                                break;
                            }
                        }
                        
                        // Add to pending list if not already there and space available
                        if (!alreadyHave && !alreadyAdded && n->grainsToAdd < 9) {
                            n->addGrainsHere[n->grainsToAdd] = grainToAdd;
                            // Find and store the orientation from the neighbor
                            for (int rg2 = 0; rg2 < nbr->grainsHere; rg2++) {
                                if (nbr->activeGrains[rg2] == grainToAdd) {
                                    n->addGrainsOrientations[n->grainsToAdd] = nbr->orientations[rg2];
                                    break;
                                }
                            }
                            n->grainsToAdd++;
                        }
                    }
                }
            }
        }
        if (n->temp < (mConfig.meltTemp)) {
            grainExists = 0;
            if (n->grainsHere > 0 ) grainExists = 1;
            if (n->grainsToAdd > 0) grainExists = 1;
            if (!grainExists && prob_dist(gen) < (1 - exp(pow(mConfig.dx,3) * mConfig.dt * n->calcNucRate(mConfig) * -1)) && n->phase < 0.1) { 
                 n->homoNucleateHere = true;
             }
            if (!grainExists &&  n->temp < mConfig.meltTemp-(mConfig.hetNucUnderCooling)&& prob_dist(gen) < n->particleComp*mConfig.dt && !n->hasHetNucleated && n->phase < 0.1) { 
                 n->hetNucleateHere = true;
                 std::cout << "Node " << n->id << " eligible for heterogeneous nucleation. Undercooling: " << (((1000*(-0.0212 *n->temp + 61.3952)/(mConfig.molarVolume))*1e-6) ) << " J/m^3, Particle Comp: " << n->particleComp << std::endl;
             }
         }
        
        // Update temperature after checking for nucleation
    n->temp = tGrad[ptr];
    }
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        globAvailPart += n->particleComp*(1-n->phase);
        if (!n) continue;
        // Add new grains from neighbors
        for (int ag = 0; ag < n->grainsToAdd; ag++) {
            if (n->grainsHere < 9) {
                int newGrainID = n->addGrainsHere[ag];
                n->activeGrains[n->grainsHere] = newGrainID;
                // Use the pre-stored orientation from addGrainsOrientations
                n->orientations[n->grainsHere] = n->addGrainsOrientations[ag];
                n->grainsHere++;
            }
        }
        n->grainsToAdd = 0;
        n->addGrainsHere = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
        n->addGrainsOrientations = {};
        if (n->hetNucleateHere) {
            // Handle heterogeneous nucleation
            n->hetNucleateHere = false;
            n->nucleatedHere = true;
            n->hasHetNucleated = true;  // Mark node as having undergone het nucleation (permanent)
            n->hetNucleationSite = true;  // Track location for visualization (red)
            globalField.addGrain(n);
            
            float sumChanged = 0;
            int distance = 1;
            std::cout << "Heterogeneous Nucleation at Node " << n->id << " (x=" << n->xPos << ", y=" << n->yPos << ")\n";
          

          
        }
        if (n->homoNucleateHere) {  
            n->nucleatedHere = true;
            n->homoNucleationSite = true;  // Track location for visualization (pink)
            // Handle homogeneous nucleation
            n->homoNucleateHere = false;
            globalField.addGrain(n);
            
            // Calculate total liquid fraction across all nodes
            
            
            
            std::cout << "Homogeneous Nucleation at Node " << n->id << " (x=" << n->xPos << ", y=" << n->yPos << ")\n";
            
        }

        // Handle homogeneous nucleation
    }
    auto t3_end = std::chrono::steady_clock::now();
    if (enableProfiling) {
        totalPropagateGrainTime += std::chrono::duration_cast<std::chrono::microseconds>(t3_end - t3_start).count();
    }

    // Report performance every 1000 updates
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
            
            // Reset counters
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


float node::calcNucRate(config modelConf) {

    double k  = 1.380649e-23;        // Boltzmann constant (J/K)
    double h  = 6.62607015e-34;      // Planck constant (J·s)
    double NA = 6.02214076e23;       // Avogadro (1/mol)

    // Driving force (J/m^3) from J/mol to J/m^3
    double dForceMo = (modelConf.drivingForceSlopek * temp
                     + modelConf.drivingForceIntercept)*5000/modelConf.molarVolume;
    

    // === 3D modification: use cell *volume*, not area ===
    double cellVolume = pow(modelConf.dx,3);

    // atomic number density (atoms / m^3)
    double atomicDensity = NA / modelConf.molarVolume;

    // atoms in the control volume
    double numAtoms = atomicDensity * cellVolume;

    // Interfacial energy γ (J/m^2)
    double gamma = 0.05;

    // === 3D spherical nucleus ΔG* ===
    double dGstar = (16.0 * 3.141592 * pow(gamma, 3.0)) /
                    (3.0 * dForceMo * dForceMo);

    // thermal factor exp(−Q / RT)
    double diffTerm = exp(-modelConf.diffusionActivationEnergy /
                          (8.314462618 * temp));

    // nucleation prefactor A = (kT / h) * N
    double prefactor = (k * temp / h) * numAtoms;

    // CNT rate: I = A * exp(−ΔG* / kT)
    double rate =prefactor * diffTerm * exp(-dGstar / (k * temp));
    //std::cout << "Nucleation Rate at Temp " << temp << " K: " << rate << " nuclei/m^3/s\n";
    return (float)rate;
}