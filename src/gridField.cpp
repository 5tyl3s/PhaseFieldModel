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

    buildGrid();
    top.grainPhases = {};
    top.phase = 0.00;
    top.exists = 0;
    bottom.grainPhases = {};
    bottom.phase = 1.00;
    bottom.exists = 0;
    // mark slots as empty with -1 (valid grain IDs start at 0)
    top.activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    bottom.activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    numGrains = 0;

    for (int ik = 0; ik < GRID_ROWS; ik++) {
        for (int jk = 0; jk < GRID_COLS; jk++) {
            grid[ik][jk].phase = 0.0;
            grid[ik][jk].particleComp = modelConfig.particleVolFraction;
            grid[ik][jk].exists = 1;
            grid[ik][jk].grainsHere = 0;
            grid[ik][jk].grainPhases = {0,0,0,0,0,0,0,0,0};
            grid[ik][jk].activeGrains = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
            grid[ik][jk].id = ik * GRID_COLS + jk;
            grid[ik][jk].heightPos = ik;
            // Initialize temperature: startTemp + vertical temperature gradient
            grid[ik][jk].temp = modelConfig.startTemp + (ik * modelConfig.tGrad * modelConfig.dx);
        };
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
    nucleus->grainsHere = nucleus->grainsHere + 1;
    nucleus->activeGrains[nucleus->grainsHere-1] = numGrains-1;
    nucleus->grainPhases[nucleus->grainsHere-1] = 1;
    nucleus->phase = 1;
    nucleus->orientations[nucleus->grainsHere-1] = tempRots;

    nucleus->neighbors[0]->grainsHere = nucleus->neighbors[0]->grainsHere + 1;
    nucleus->neighbors[0]->activeGrains[nucleus->neighbors[0]->grainsHere-1] = numGrains-1;
    nucleus->neighbors[0]->orientations[nucleus->neighbors[0]->grainsHere-1] = tempRots;
    nucleus->neighbors[1]->grainsHere = nucleus->neighbors[1]->grainsHere + 1;
    nucleus->neighbors[1]->activeGrains[nucleus->neighbors[1]->grainsHere-1] = numGrains-1;
    nucleus->neighbors[1]->orientations[nucleus->neighbors[1]->grainsHere-1] = tempRots;




    if (nucleus->neighbors[2]->exists != 0) {
        nucleus->neighbors[2]->grainsHere = nucleus->neighbors[2]->grainsHere + 1;
        nucleus->neighbors[2]->activeGrains[nucleus->neighbors[2]->grainsHere-1] = numGrains-1;
        nucleus->neighbors[2]->orientations[nucleus->neighbors[2]->grainsHere-1] = tempRots;
    }
    if (nucleus->neighbors[3]->exists != 0) {
        nucleus->neighbors[3]->grainsHere = nucleus->neighbors[3]->grainsHere + 1;
        nucleus->neighbors[3]->activeGrains[nucleus->neighbors[3]->grainsHere-1] = numGrains-1;
        nucleus->neighbors[3]->orientations[nucleus->neighbors[3]->grainsHere-1] = tempRots;
    }



}





void gridField::update(
    std::array<double, TOTAL_NODES> &phaseDiffEn,
    std::array<std::array<double,9>, TOTAL_NODES> &grainDiffEn,
    std::array<double, TOTAL_NODES> &tempPartComp, // now holds chemical potential μ at each node
    std::array<double, TOTAL_NODES> &tGrad
) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    bool flip;
    bool grainExists = 0;

    // use compile-time TOTAL_NODES from header
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        // make sure pointer is valid (defensive)
        node* n = allNodes[ptr];
        if (!n) continue;
        
        // update phase
        n->phase = n->phase -(mConfig.dt / pow(mConfig.dx,2)) * phaseDiffEn[ptr]*5e-20;

        // update grain phases
        for (int g = 0; g < 9; g++) {
            n->grainPhases[g] = n->grainPhases[g] - (mConfig.dt / pow(mConfig.dx,2) * grainDiffEn[ptr][g]*1e-19);
        }
    }

    // THERMODYNAMIC particle update: dC/dt = mobility * laplacian(μ)
    // Particles flow DOWN the chemical potential gradient to minimize free energy
    // μ is computed in calcParticleCompDiff() and stored in tempPartComp array
    const double particleMobilityLiquid = 5e-19;   // high mobility in liquid
    const double particleMobilitySolid = 1e-28;    // lower mobility in solid
    double dx = mConfig.dx;
    double denom = (dx * dx);
    
    // Store composition changes
    std::array<double, TOTAL_NODES> particleCompChange;
    particleCompChange.fill(0.0);
    
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        if (!n) continue;

        // tempPartComp[ptr] contains the chemical potential μ at this node
        double mu = tempPartComp[ptr];
        
        // Compute 4-point Laplacian of chemical potential
        double sumNeighMu = 0.0;
        int nExist = 0;
        for (int nb = 0; nb < 4; nb++) {  // ONLY 4-point: left, right, up, down
            node* nbr = n->neighbors[nb];
            if (!nbr) continue;
            if (nbr->exists == 0) continue; // no-flux boundary
            sumNeighMu += tempPartComp[nbr->id];
            nExist++;
        }

        double lapMu = 0.0;
        if (nExist > 0) {
            lapMu = (sumNeighMu - nExist * mu);  // (sum of 4 neighbors - 4*center)
        }

        // Choose mobility based on phase (liquid vs solid)
        double liquidFrac = (1.0 - n->phase);  // 1 in liquid, 0 in solid
        double mobility = particleMobilityLiquid * liquidFrac + particleMobilitySolid * (1.0 - liquidFrac);
        
        // Flux: J = mobility * laplacian(μ)
        // This drives particles toward lower μ regions (energy minima)
        double flux = mobility * lapMu / denom;
        
        // Store change
        particleCompChange[ptr] = -1*mConfig.dt * flux;
    }
    
    // Apply changes TOGETHER to preserve global sum
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        if (!n) continue;
        n->particleComp += particleCompChange[ptr];
        
        // Clamp to [0, 1] to prevent unphysical values
        if (n->particleComp < 0.0) n->particleComp = 0.0;
        if (n->particleComp > 1.0) n->particleComp = 1.0;
    }

    // second pass: clamp values and propagate active grains / nucleation
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        if (!n) continue;

        if (n->phase > 1.0) n->phase = 1.0;
        if (n->phase < 0.0) n->phase = 0.0;

        for (int g = 0; g < n->grainsHere; g++) {
            int activeGrainID = n->activeGrains[g];
            if (activeGrainID < 0) continue;

            // Clamp the local grain slot value (index = g)
            if (n->grainPhases[g] >= 1.0) {
                n->grainPhases[g] = 1.0;
            }
            if (n->grainPhases[g] < 0.0) {
                n->grainPhases[g] = 0.0;
            }

            if (n->grainPhases[g] > 1e-16) {
                // try to add this grain (by global id) to neighbors if not already present
                for (int nb = 0; nb < 8; nb++) {
                    node* nbr = n->neighbors[nb];
                     if (!nbr || nbr->exists == 0) continue;
                     bool found = false;
                     for (int rg = 0; rg < nbr->grainsHere; rg++) {
                         if (nbr->activeGrains[rg] == activeGrainID) {
                             found = true;
                             break;
                         }
                     }
                     if (!found) {
                         int insertPos = nbr->grainsHere;
                         if (insertPos < 9) {
                             nbr->grainsHere = insertPos + 1;
                             nbr->activeGrains[insertPos] = activeGrainID;
                             nbr->orientations[insertPos] = n->orientations[g];
                         }
                     }
                 }
             }
         }

        if (n->temp < (mConfig.meltTemp)) {
            grainExists = 0;
            if (n->grainsHere > 0 ) grainExists = 1;
            for (int nb = 0; nb < 8 && !grainExists; nb++) {
                node* nbr = n->neighbors[nb];
                if (!nbr) continue;
                if (nbr->grainsHere > 0) grainExists = 1;
            }
             
            if (!grainExists && prob_dist(gen) < (1 - exp(pow(mConfig.dx,3) * mConfig.dt * n->calcNucRate(mConfig) * -1))) { 
                 addGrain(n);
             }
         }
        
        // Update temperature after checking for nucleation
        n->temp = tGrad[ptr];
             
         
     }

}


float node::calcNucRate(config modelConf) {

    double k  = 1.380649e-23;        // Boltzmann constant (J/K)
    double h  = 6.62607015e-34;      // Planck constant (J·s)
    double NA = 6.02214076e23;       // Avogadro (1/mol)

    // Driving force (J/m^3)
    double dForceMo = modelConf.drivingForceSlopek * temp
                    + modelConf.drivingForceIntercept;
    dForceMo = 1000.0 * dForceMo / modelConf.molarVolume;

    // === 3D modification: use cell *volume*, not area ===
    double cellVolume = pow(modelConf.dx,3);

    // atomic number density (atoms / m^3)
    double atomicDensity = NA / modelConf.molarVolume;

    // atoms in the control volume
    double numAtoms = atomicDensity * cellVolume;

    // Interfacial energy γ (J/m^2)
    double gamma = modelConf.liqSolIntE;

    // === 3D spherical nucleus ΔG* ===
    double dGstar = (16.0 * 3.141592 * pow(gamma, 3.0)) /
                    (3.0 * dForceMo * dForceMo);

    // thermal factor exp(−Q / RT)
    double diffTerm = exp(-modelConf.diffusionActivationEnergy /
                          (8.314462618 * temp));

    // nucleation prefactor A = (kT / h) * N
    double prefactor = (k * temp / h) * numAtoms;

    // CNT rate: I = A * exp(−ΔG* / kT)
    double rate = prefactor * diffTerm * exp(-dGstar / (k * temp));
    //std::cout << "Nucleation Rate at Temp " << temp << " K: " << rate << " nuclei/m^3/s\n";
    return (float)rate;
}