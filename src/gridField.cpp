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
        gp.neighbors[0] = (j > 0) ? &grid[0][j - 1] : &grid[0][GRID_COLS - 1]; // Left
        gp.neighbors[1] = (j < GRID_COLS - 1) ? &grid[0][j + 1] : &grid[0][0]; // Right
        gp.neighbors[2] = &top;
        gp.neighbors[3] = &grid[1][j]; // Down
        allNodes[j] = &gp;
        gp.heightPos = 0;
    }

    // interior rows
    for (int i = 1; i < GRID_ROWS - 1; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            node& gp = grid[i][j];
            gp.neighbors[0] = (j > 0) ? &grid[i][j - 1] : &grid[i][GRID_COLS - 1]; // Left
            gp.neighbors[1] = (j < GRID_COLS - 1) ? &grid[i][j + 1] : &grid[i][0]; // Right
            gp.neighbors[2] = &grid[i - 1][j]; // Up
            gp.neighbors[3] = &grid[i + 1][j]; // Down
            allNodes[i * GRID_COLS + j] = &gp;
            gp.heightPos = i;
        }
    }

    // bottom row
    for (int j = 0; j < GRID_COLS; j++) {
        node& gp = grid[GRID_ROWS - 1][j];
        gp.neighbors[0] = (j > 0) ? &grid[GRID_ROWS - 1][j - 1] : &grid[GRID_ROWS - 1][GRID_COLS - 1]; // Left
        gp.neighbors[1] = (j < GRID_COLS - 1) ? &grid[GRID_ROWS - 1][j + 1] : &grid[GRID_ROWS - 1][0]; // Right
        gp.neighbors[2] = &grid[GRID_ROWS - 2][j];
        gp.neighbors[3] = &bottom; // Down
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
        

        n->temp = tGrad[ptr];

        // update phase
        n->phase = n->phase - (mConfig.dt / mConfig.cellArea) * phaseDiffEn[ptr];

        // update grain phases
        for (int g = 0; g < 9; g++) {
            n->grainPhases[g] = n->grainPhases[g] - (1e15 * mConfig.dt / mConfig.cellArea * grainDiffEn[ptr][g]);
        }
    }

    // conservative particle composition update using chemical potential array (tempPartComp)
    // discrete divergence of flux: dC/dt = mobility * laplacian(mu)
    const double particleMobility = 1e-6; // tune or add to config
    double dx = mConfig.dx;
    double denom = (dx * dx);
    for (int ptr = 0; ptr < TOTAL_NODES; ptr++) {
        node* n = allNodes[ptr];
        if (!n) continue;

        double mu = tempPartComp[ptr];

        double sumNeighMu = 0.0;
        int nExist = 0;
        for (int nb = 0; nb < 4; nb++) {
            node* nbr = n->neighbors[nb];
            if (!nbr) continue;
            if (nbr->exists == 0) continue; // enforce no-flux by skipping
            sumNeighMu += tempPartComp[nbr->id];
            nExist++;
        }

        double lapMu = 0.0;
        if (nExist > 0) {
            lapMu = (sumNeighMu - nExist * mu);
        }

        // Compute change in composition (conservative)
        double dC =  (mConfig.dt / denom) * lapMu;

        // Optionally gate diffusion by liquid fraction (keep particles mostly mobile in liquid)
        double mobilityGate = (1.0 - n->phase); // 1 in liquid, 0 in fully solid
        n->particleComp += dC * mobilityGate;

        // clamp to physical bounds
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

            if (n->grainPhases[g] > 1e-9) {
                // try to add this grain (by global id) to neighbors if not already present
                for (int nb = 0; nb < 4; nb++) {
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
            if (n->grainsHere > 0 || n->phase > 0.000000001 || n->neighbors[0]->grainsHere > 0 || n->neighbors[1]->grainsHere > 0 || n->neighbors[2]->grainsHere > 0 || n->neighbors[3]->grainsHere > 0) {
                grainExists = 1;
            }
            
            if (!grainExists && prob_dist(gen) < (1 - exp(mConfig.dx * mConfig.dx * mConfig.dt * n->calcNucRate(mConfig) * -1))) { 
                addGrain(n);
            }
        }
            
        
    }

}


float node::calcNucRate(config modelConf) {
    double dForceMo = modelConf.drivingForceSlopek * temp + modelConf.drivingForceIntercept;
    // kJ/mol
    dForceMo = 1000*dForceMo / modelConf.molarVolume;
    // J/m^3
    //boltzman constant
    //std::cout << "Driving Force: " << dForceMo << std::endl;

    double k= 1.380649e-23; // m^2kg/(s^2K)
    // planck constant
    double h = 6.62607015e-34; // m^2kg/s
    //Number atoms in cell
    double avogadro = 6.02214076e23; // #atoms/mole
    // (m^3/mol)^(2/3) * (#atoms/mol)^(1/3)
    double numAtoms = modelConf.cellArea*avogadro/(pow(modelConf.molarVolume, 2.0/3.0) * pow(avogadro, 1.0/3.0));
    //std::cout << "Number of Atoms in Cell: " << numAtoms << std::endl;
    //std::cout << "LiqSolIntE: " << modelConf.liqSolIntE << std::endl;
    double freeEnergyNucleusFormation = (16*3.14159*pow(modelConf.liqSolIntE,3)/(3*(dForceMo*dForceMo)));
    // 16*pi*<liqSolIntE^3>/(3*(dForceMo^2)) = (J/m^2)^3 / (J/m^3)^2 = J
    //
    //std::cout << "Free Energy of Nucleus Formation at Temp " << temp << " is " << freeEnergyNucleusFormation << std::endl;
    long double inner = freeEnergyNucleusFormation/(k*temp);
    // J / (m^2kg/(s^2K) * K) = J/(m^2kg/(s^2)) = (kg*m^2/s^2)/(m^2kg/s^2) = unitless

    //std::cout << "Inner Exponential Term: " << std::setprecision(10) << inner << std::endl;
    long double nucleationRate = (k*temp*numAtoms/h)*exp(-1*modelConf.diffusionActivationEnergy/(8.314*temp))*exp(-inner);
    //std::cout << "First Term: "<< (k*temp*numAtoms/h) << std::endl << "Second Term: " << exp(-1*modelConf.diffusionActivationEnergy/(8.314*temp)) << std::endl << "ThirdTerm: " << exp(-1*inner) << std::endl;
   // std::cout << "Nucleation Rate at Temp " << temp << " is " << nucleationRate << std::endl << std::endl;

    return nucleationRate;
}