#include "gridField.hpp"
#include <iostream>
#include <random>
#include <cmath>
// This line actually creates the global object in memory
gridField globalField;
std::array<std::array<double,1000>,1000> tempGrid;

void gridField::buildGrid() {


    
    
    for (int j = 0; j< 1000;j++) {
        node& gp = grid[0][j];
        gp.neighbors[0] = (j>0) ? &grid[0][j-1]:&grid[0][1000-1]; //Left
        gp.neighbors[1] = (j<1000-1) ? &grid[0][j+1]:&grid[0][0]; //Right
        gp.neighbors[2] =&top;
        gp.neighbors[3] = &grid[1][j]; //Down
        allNodes[j] = &gp;
        gp.heightPos =0;
        
    }
    for (int i = 1; i< 1000-1;i++) {
        for (int j = 0; j < 1000;j++) {
            node& gp = grid[i][j];
            gp.neighbors[0] = (j>0) ? &grid[i][j-1]:&grid[i][1000-1]; //Left
            gp.neighbors[1] = (j<1000-1) ? &grid[i][j+1]:&grid[i][0]; //Right
            gp.neighbors[2] = (i>0) ? &grid[i-1][j]: nullptr;//up
            gp.neighbors[3] = (i<1000-1) ? &grid[i+1][j]: nullptr; //Down
            allNodes[i*1000 + j] = &gp;
            gp.heightPos = i;
        };
    };
    for (int j = 0; j< 1000;j++) {
        node& gp = grid[1000-1][j];
        gp.neighbors[0] = (j>0) ? &grid[1000-1][j-1]:&grid[1000-1][1000-1]; //Left
        gp.neighbors[1] = (j<1000-1) ? &grid[1000-1][j+1]:&grid[1000-1][0]; //Right
        gp.neighbors[2] = &grid[1000-2][j];
        gp.neighbors[3] = &bottom; //Down
        allNodes[(1000-1)*1000 + j] = &gp;
        gp.heightPos = (1000-1);
    }

};


void gridField::init(config modelConfig) {

    


    buildGrid();
    top.grainPhases = {};
    top.phase = 0.00;
    top.exists = 0;
    bottom.grainPhases = {};
    bottom.phase = 1.00;
    bottom.exists = 0;
    numGrains = 0;

    for (int ik = 0; ik < 1000; ik++) {
        for (int jk = 0; jk < 1000; jk++) {
            
            grid[ik][jk].phase = 0.0;
            grid[ik][jk].particleComp = modelConfig.particleVolFraction; 
            grid[ik][jk].exists = 1;
            grid[ik][jk].grainsHere = 0;
            grid[ik][jk].grainPhases = {0,0,0,0,0,0,0,0,0};
            grid[ik][jk].activeGrains = {0,0,0,0,0,0,0,0,0};
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
    std::array<double,1000000> &phaseDiffEn,
    std::array<std::array<double,9>,1000000> &grainDiffEn,
    std::array<double,1000000> &tempPartComp,
    std::array<double,1000000> &tGrad
) {
    //std::cout << "\n\n\nStart Update:\n"; 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    bool flip;
    bool grainExists = 0;


    for (int ptr = 0; ptr < 1000000; ptr++) {
        allNodes[ptr]->temp = tGrad[ptr];
        allNodes[ptr]->phase = allNodes[ptr]->phase - mConfig.dt /mConfig.cellArea * phaseDiffEn[ptr];
        allNodes[ptr]->particleComp = allNodes[ptr]->particleComp - ((1/mConfig.cellArea)*mConfig.dt*(1-allNodes[ptr]->phase)*tempPartComp[ptr]);
        for (int g = 0; g < 9; g++) {
            allNodes[ptr]->grainPhases[g] = allNodes[ptr]->grainPhases[g] - (1e15*mConfig.dt/mConfig.cellArea*grainDiffEn[ptr][g]);
        }
        

    }


    
  
    //std::cout << "Starting Grain Activation Update" << std::endl;   

 


    for (int ptr = 0; ptr < 1000000; ptr++) {

            if (allNodes[ptr]->phase > 1) {
                allNodes[ptr]->phase = 1;
            }
            if (allNodes[ptr]->phase < 0) {
                allNodes[ptr]->phase = 0;
            }
            for (int g = 0; g < allNodes[ptr]->grainsHere; g++) {
                //std::cout << "Grain Activation Update: (" << i << "," << j << ") Grain " << grid[i][j].activeGrains[g] << " Phase: " << grid[i][j].grainPhases[grid[i][j].activeGrains[g]] << std::endl;
                if (allNodes[ptr]->grainPhases[allNodes[ptr]->activeGrains[g]] >= 1) {
                    allNodes[ptr]->grainPhases[allNodes[ptr]->activeGrains[g]] = 1;
                }
                if (allNodes[ptr]->grainPhases[g] < 0) {
                        allNodes[ptr]->grainPhases[g] = 0;
                }
                if (allNodes[ptr]->grainPhases[g] > 0.000000001) {
                    if (allNodes[ptr]->neighbors[0]->exists != 0) {
                        flip = 0;
                        for (int rg = 0; rg < allNodes[ptr]->neighbors[0]->grainsHere; rg++) {
                            if (allNodes[ptr]->neighbors[0]->activeGrains[rg] == allNodes[ptr]->activeGrains[g]) {
                                flip = 1;   
                            
                            }
                        }
                        if (flip == 0) {

                            allNodes[ptr]->neighbors[0]->grainsHere = allNodes[ptr]->neighbors[0]->grainsHere + 1;
                            allNodes[ptr]->neighbors[0]->activeGrains[allNodes[ptr]->neighbors[0]->grainsHere-1] = allNodes[ptr]->activeGrains[g];
                            allNodes[ptr]->neighbors[0]->orientations[allNodes[ptr]->neighbors[0]->grainsHere-1] = allNodes[ptr]->orientations[g];
                        }
                    }

                    if (allNodes[ptr]->neighbors[1]->exists != 0) {
                        flip = 0;
                        for (int rg = 0; rg < allNodes[ptr]->neighbors[1]->grainsHere; rg++) {
                            if (allNodes[ptr]->neighbors[1]->activeGrains[rg] == allNodes[ptr]->activeGrains[g]) {
                                flip = 1;   
                            
                            }
                        }
                        if (flip == 0) {

                            allNodes[ptr]->neighbors[1]->grainsHere = allNodes[ptr]->neighbors[1]->grainsHere + 1;
                            allNodes[ptr]->neighbors[1]->activeGrains[allNodes[ptr]->neighbors[1]->grainsHere-1] = allNodes[ptr]->activeGrains[g];
                            allNodes[ptr]->neighbors[1]->orientations[allNodes[ptr]->neighbors[1]->grainsHere-1] = allNodes[ptr]->orientations[g];
                        }
                    }
                    if (allNodes[ptr]->neighbors[2]->exists != 0) {
                        flip = 0;
                        for (int rg = 0; rg < allNodes[ptr]->neighbors[2]->grainsHere; rg++) {
                            if (allNodes[ptr]->neighbors[2]->activeGrains[rg] == allNodes[ptr]->activeGrains[g]) {
                                flip = 1;   
                            
                            }
                        }
                        if (flip == 0) {

                            allNodes[ptr]->neighbors[2]->grainsHere = allNodes[ptr]->neighbors[2]->grainsHere + 1;
                            allNodes[ptr]->neighbors[2]->activeGrains[allNodes[ptr]->neighbors[2]->grainsHere-1] = allNodes[ptr]->activeGrains[g];
                            allNodes[ptr]->neighbors[2]->orientations[allNodes[ptr]->neighbors[2]->grainsHere-1] = allNodes[ptr]->orientations[g];
                        }
                    }
                    if (allNodes[ptr]->neighbors[3]->exists != 0) {
                        flip = 0;
                        for (int rg = 0; rg < allNodes[ptr]->neighbors[3]->grainsHere; rg++) {
                            if (allNodes[ptr]->neighbors[3]->activeGrains[rg] == allNodes[ptr]->activeGrains[g]) {
                                flip = 1;   
                            
                            }
                        }
                        if (flip == 0) {

                            allNodes[ptr]->neighbors[3]->grainsHere = allNodes[ptr]->neighbors[3]->grainsHere + 1;
                            allNodes[ptr]->neighbors[3]->activeGrains[allNodes[ptr]->neighbors[3]->grainsHere-1] = allNodes[ptr]->activeGrains[g];
                            allNodes[ptr]->neighbors[3]->orientations[allNodes[ptr]->neighbors[3]->grainsHere-1] = allNodes[ptr]->orientations[g];
                        }
                }
            }
        



            
            if (allNodes[ptr]->temp < (mConfig.meltTemp)) {
                grainExists = 0;
                if (allNodes[ptr]->grainsHere > 0 || allNodes[ptr]->phase > 0.000000001 || allNodes[ptr]->neighbors[0]->grainsHere > 0 || allNodes[ptr]->neighbors[1]->grainsHere > 0 || allNodes[ptr]->neighbors[2]->grainsHere > 0 || allNodes[ptr]->neighbors[3]->grainsHere > 0) {
                    grainExists = 1;
                }


                if (!grainExists && prob_dist(gen) <(1-exp(mConfig.dx*mConfig.dx*mConfig.dt*allNodes[ptr]->calcNucRate(mConfig)*-1))) { 

                    addGrain(allNodes[ptr]);
                }
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
    long double nucleationRate = (k*temp*numAtoms/h)*exp(-1*modelConf.diffusionActivationEnergy/(8.314*temp))*exp(inner);
    //std::cout << "First Term: "<< (k*temp*numAtoms/h) << std::endl << "Second Term: " << exp(-1*modelConf.diffusionActivationEnergy/(8.314*temp)) << std::endl << "ThirdTerm: " << exp(-1*inner) << std::endl;
   // std::cout << "Nucleation Rate at Temp " << temp << " is " << nucleationRate << std::endl << std::endl;

    return nucleationRate;
}