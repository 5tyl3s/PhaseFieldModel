#include <iostream>
#include "modelStartup.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <cmath>
#include "json.hpp" 
using json = nlohmann::json;

config inputConfig() {
    config newConfig;
    std::string configPath = "../config/modelConfig.json";
    std::ifstream configIn(configPath);
    if (!configIn.is_open()) {
        std::cout << "Model Config File Not Found. Generating default config at: " << configPath << std::endl;
        std::ofstream configOut(configPath);
        if (!configOut) {
            std::cerr << "Failed to create config file!" << std::endl;
            newConfig.success =  0;
            return newConfig;
        }
        configOut << "{\n";
        configOut << "  \"dx_um\": 0.5,\n";
        configOut << "  \"meltTemp\": 2896,\n";
        configOut << "  \"HeatCapacityJKgK\": 421,\n";
        configOut << "  \"Densitykgm3\": 10280,\n";
        configOut << "  \"dt_s\": 0.001,\n";
        configOut << "  \"startTempK\": 3000,\n";
        configOut << "  \"solidConductivityWMK\": 138,\n";
        configOut << "  \"cellHeightum\": 100,\n";
        configOut << "  \"cellWidthum\": 100,\n";
        configOut << "  \"particleVolumeFraction\": 0.01,\n";
        configOut << "  \"LiqSolIntEnergyJm2\": 0.3,\n";
        configOut << "  \"LiqSolIntWidthum\": 0.3,\n";
        configOut << "  \"GrainIntEnergyJm2\": 0.3,\n";
        configOut << "  \"GrainIntWidthum\": 0.3,\n";
        configOut << "  \"PhaseBarrierHeightCoefficient\": 0.25,\n";
        configOut << "  \"BasePlateTempK\": 500.00,\n";
        configOut << "  \"GrainBarrierHeightCoefficient\": 0.125\n";

        configOut << "}\n";
        configOut.close();
        std::cout << "Default Config Created. Please Edit and Rerun Program";
        newConfig.success = 0;


    }
    std::string line;
    json j;
    configIn >> j;

    newConfig.dx = j["dx_um"];
    newConfig.dx = newConfig.dx*1e-6;
    newConfig.dt = j["dt_s"];
    newConfig.startTemp = j["startTempK"];
    newConfig.kSolid = j["solidConductivityWMK"];
    newConfig.cellHeight = j["cellHeightum"];
    newConfig.cellHeight = newConfig.cellHeight*1e-6;
    newConfig.cellWidth = j["cellWidthum"];
    newConfig.cellWidth = newConfig.cellWidth*1e-6;
    newConfig.cellArea = newConfig.cellHeight * newConfig.cellWidth;
    newConfig.particleVolFraction = j["particleVolumeFraction"];
    newConfig.liqSolIntE = j["LiqSolIntEnergyJm2"];
    newConfig.liqSolIntWidth = j["LiqSolIntWidthum"];
    newConfig.phaseCoefficient = 0.75*newConfig.liqSolIntWidth * newConfig.liqSolIntE;
    newConfig.barrierHeightGrain = j["GrainBarrierHeightCoefficient"];
    newConfig.barrierHeightPhase = j["PhaseBarrierHeightCoefficient"];
    newConfig.basePlateTemp = j["BasePlateTempK"];
    newConfig.meltTemp = j["meltTemp"];
    newConfig.heatCapacity = j["HeatCapacityJKgK"];
    newConfig.density = j["Densitykgm3"];
    newConfig.kLiquid = newConfig.kSolid * 0.45; 
    newConfig.phasePreCo = 0.75*newConfig.liqSolIntE/(newConfig.barrierHeightPhase*newConfig.liqSolIntWidth);
    newConfig.grainPreCo = 0.75*newConfig.grainIntE/(newConfig.barrierHeightGrain*newConfig.grainIntWidth);
    newConfig.phaseGradCo = 0.75*newConfig.liqSolIntE*newConfig.liqSolIntWidth;
    newConfig.particleSlowingCoefficient = 0;
    newConfig.success = 1;
    
 



    std::cout << "The Step Size is " << newConfig.dx << std::endl;
    int height = static_cast<int>(std::ceil(newConfig.cellHeight / newConfig.dx));



    int width = static_cast<int>(std::ceil(newConfig.cellWidth / newConfig.dx));
    std::vector<int> steps = {height,width};
    newConfig.steps = steps;


    int totalSteps = steps[0] * steps[1];
    newConfig.totalSteps = totalSteps;

   
    
    

    newConfig.success = 1;
    return newConfig;
};



void gridField::buildGrid(int widthSteps, int heightSteps) {
    std::cout << widthSteps << " wide " << heightSteps << " Tall" << std::endl;
    grid.resize(heightSteps,std::vector<node>(widthSteps));


 
  
    

    for (int j = 0; j< widthSteps;j++) {
        node& gp = grid[0][j];
        gp.neighbors[0] = (j>0) ? &grid[0][j-1]:&grid[0][widthSteps-1]; //Left
        gp.neighbors[1] = (j<widthSteps-1) ? &grid[0][j+1]:&grid[0][0]; //Right
        gp.neighbors[2] =&top;
        gp.neighbors[3] = &grid[1][j]; //Down
        
    }

    for (int i = 1; i< heightSteps;i++) {
        for (int j = 0; j < widthSteps;j++) {
            node& gp = grid[i][j];
            gp.neighbors[0] = (j>0) ? &grid[i][j-1]:&grid[i][widthSteps-1]; //Left
            gp.neighbors[1] = (j<widthSteps-1) ? &grid[i][j+1]:&grid[i][0]; //Right
            gp.neighbors[2] = (i>0) ? &grid[i-1][j]: nullptr;//up
            gp.neighbors[3] = (i<heightSteps-1) ? &grid[i+1][j]: nullptr; //Down
        };
    };

    for (int j = 0; j< widthSteps;j++) {
        node& gp = grid[heightSteps-1][j];
        gp.neighbors[0] = (j>0) ? &grid[heightSteps-1][j-1]:&grid[heightSteps-1][widthSteps-1]; //Left
        gp.neighbors[1] = (j<widthSteps-1) ? &grid[heightSteps-1][j+1]:&grid[heightSteps-1][0]; //Right
        gp.neighbors[2] = &grid[heightSteps-2][j];
        gp.neighbors[3] = &bottom; //Down
        
    }





};


void gridField::init(config modelConfig) {
    

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0,1);
    buildGrid(modelConfig.steps[0],modelConfig.steps[1]);
    top.grainPhases = {};
    top.phase = 0.00;
    bottom.grainPhases = {};
    bottom.phase = 1.00;
    numGrains = 0;

    for (int ik = 0; ik < modelConfig.steps[0]; ik++) {
        for (int jk = 0; jk < modelConfig.steps[1]; jk++) {
            grid[ik][jk].temp = modelConfig.startTemp +dist(gen);
            grid[ik][jk].phase = 0.0;
            grid[ik][jk].particleComp = modelConfig.particleVolFraction; 
        };
    }
}


void gridField::addGrain(std::vector<int> nucleus,config modelConf) {
    
    numGrains = numGrains + 1;
    std::cout << "There are now "<<numGrains << "Grains\n";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(45,45);

    for (int i = 0; i < modelConf.steps[0]; i++) {
        for (int j = 0; j < modelConf.steps[1]; j++) {
            grid[i][j].grainPhases.push_back(0.0);
        }
    }
    


    top.grainPhases.push_back(0.0);
    bottom.grainPhases.push_back(0.0);
    
    double one = 1.000;
    //std::cout << "hi";

    grid[nucleus[0]][nucleus[1]].grainPhases[numGrains-1] = one;
    //std::cout << "hi";
    eulerAngles tempRots = {dist(gen),dist(gen),dist(gen)};
    ////std::cout << tempRots.theta1;
    orientations.push_back(tempRots);
    //std::cout << "Hi";

}
void gridField::update(std::vector<std::vector<double>> phaseDiffEn, std::vector<std::vector<double>> tempGrad, std::vector<std::vector<std::vector<double>>> grainDiffEn, config modelConf) {
    //std::cout << "\n\n\nStart Update:\n"; 
    for (int i = 0; i < modelConf.steps[0]; i++) {
        for (int j = 0; j < modelConf.steps[1]; j++) {
            //std::cout <<"OldTemp: " << grid[i][j].temp;
            //std::cout << "Thermal Diffusivity: " << (((modelConf.kLiquid+(grid[i][j].phase*(modelConf.kSolid-modelConf.kLiquid)))/(modelConf.density*modelConf.heatCapacity))) << std::endl;
            //std::cout << "Keff: " << modelConf.kLiquid+(grid[i][j].phase*(modelConf.kSolid-modelConf.kLiquid)) << std::endl;
            //std::cout << modelConf.kLiquid << std::endl;
            //std::cout << modelConf.kSolid << std::endl;
            
            //std::cout << "Phase: " << grid[i][j].phase << std::endl;
            //std::cout << "GradT: " << tempGrad[i][j] << std::endl;
            //std::cout <<" OldTemp: " << grid[i][j].temp << std::endl;
            //std::cout <<  "DT: " << (((modelConf.kLiquid+(grid[i][j].phase*(modelConf.kSolid-modelConf.kLiquid)))/(modelConf.density*modelConf.heatCapacity)))*modelConf.dt*(tempGrad[i][j]) << std::endl;
            //std::cout << "1";
            grid[i][j].temp = grid[i][j].temp + (((modelConf.kLiquid+(grid[i][j].phase*(modelConf.kSolid-modelConf.kLiquid)))/(modelConf.density*modelConf.heatCapacity)))*modelConf.dt*(tempGrad[i][j]);
            // std::cout << "2";

            //std::cout <<" NewTemp: " << grid[i][j].temp << std::endl;
            //std::cout << "Old Phase: " <<grid[i][j].phase;
            grid[i][j].phase = grid[i][j].phase + phaseDiffEn[i][j];
            //std::cout << "3";
            //std::cout << "Help";
            //std::cout << " New Phase: " <<grid[i][j].phase << std::endl;
            for (int g = 0; g < numGrains;g++) {
                //std::cout << "What";
                //std::cout << numGrains << "Nums" << std::endl;
                //std::cout << grainDiffEn[i][j][g] << std::endl;
                grid[i][j].grainPhases[g] = grid[i][j].grainPhases[g] + grainDiffEn[i][j][g];
            }

            //std::cout << "Final Phase: " << grid[i][j].phase << std::endl;
        }
    }

    for (int i = 0; i < modelConf.steps[0]; i++) {
        for (int j = 0; j < modelConf.steps[1]; j++) {
            if (grid[i][j].phase > 1) {
                grid[i][j].phase = 1;
            }
            if (grid[i][j].phase < 0) {
                grid[i][j].phase = 0;
            }
            if (grid[i][j].temp < modelConf.meltTemp*0.9) {

                if (grid[i][j].phase < 0.1) {
                    grid[i][j].phase = 1;
                    addGrain({i,j},modelConf);
                }
            }
        }
    }

}


       
      