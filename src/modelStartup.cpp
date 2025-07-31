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
        configOut << "  \"BasePlateTempK\": 500,\n";
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
    newConfig.dt = j["dt_s"];
    newConfig.startTemp = j["startTempK"];
    newConfig.kSolid = j["solidConductivityWMK"];
    newConfig.cellHeight = j["cellHeightum"];
    newConfig.cellWidth = j["cellWidthum"];
    newConfig.particleVolFraction = j["particleVolumeFraction"];
    newConfig.liqSolIntE = j["LiqSolIntEnergyJm2"];
    newConfig.liqSolIntWidth = j["LiqSolIntWidthum"];
    newConfig.phaseCoefficient = 0.75*newConfig.liqSolIntWidth * newConfig.liqSolIntE;
    newConfig.barrierHeightGrain = j["GrainBarrierHeightCoefficient"];
    newConfig.barrierHeightPhase = j["PhaseBarrierHeightCoefficient"];
    newConfig.basePlateTemp = j["BasePlateTempK"];
    newConfig.meltTemp = j["meltTemp"];

    newConfig.success = 1;
 



    std::cout << "The Step Size is " << newConfig.dx << std::endl;
    int height = static_cast<int>(std::ceil(newConfig.cellHeight / newConfig.dx));



    int width = static_cast<int>(std::ceil(newConfig.cellWidth / newConfig.dx));
    std::array<int,2> steps = {height,width};
    newConfig.steps = steps;


    int totalSteps = steps[0] * steps[1];
    newConfig.totalSteps = totalSteps;

   
    
    

    newConfig.success = 1;
    return newConfig;
};



void gridField::buildGrid(int widthSteps, int heightSteps) {
    std::cout << widthSteps << " wide " << heightSteps << " Tall" << std::endl;
    grid.resize(heightSteps,std::vector<node>(widthSteps));


 
  
    

    for (int j = 0; j< widthSteps-1;j++) {
        node& gp = grid[0][j];
        gp.neighbors[0] = (j>0) ? &grid[0][j-1]:&grid[0][widthSteps-1]; //Left
        gp.neighbors[1] = (j<widthSteps-1) ? &grid[0][j+1]:&grid[0][0]; //Right
        gp.neighbors[2] =&top;
        gp.neighbors[3] = &grid[1][j]; //Down
        
    }

    for (int i = 1; i< heightSteps-2;i++) {
        for (int j = 0; j < widthSteps-1;j++) {
            node& gp = grid[i][j];
            gp.neighbors[0] = (j>0) ? &grid[i][j-1]:&grid[i][widthSteps-1]; //Left
            gp.neighbors[1] = (j<widthSteps-1) ? &grid[i][j+1]:&grid[i][0]; //Right
            gp.neighbors[2] = (i>0) ? &grid[i-1][j]: nullptr;//up
            gp.neighbors[3] = (i<heightSteps-1) ? &grid[i+1][j]: nullptr; //Down
        };
    };

    for (int j = 0; j< widthSteps-1;j++) {
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
    std::normal_distribution<> dist(0.0,25);
    buildGrid(modelConfig.steps[0],modelConfig.steps[1]);
    top.grainPhases = {};
    top.phase = 0.00;
    bottom.grainPhases = {};
    bottom.phase = 1.00;

    for (int i = 0; i < modelConfig.steps[0]; i++) {
        for (int j = 0; j < modelConfig.steps[1]; j++) {
            grid[i][j].temp = modelConfig.startTemp +dist(gen);
            grid[i][j].phase = 0;
            grid[i][j].particleComp = modelConfig.particleVolFraction; 
        };
    }
}


void gridField::addGrain(std::array<int,2> nucleus) {
    numGrains = numGrains + 1;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(45,45);

    for (int i = 1; i < grid.size()-1; i++) {
        for (int j = 1; j < grid[1].size()-1; j++) {
            grid[i][j].grainPhases.push_back(0);
        }
    }
    top.grainPhases.push_back(0);
    bottom.grainPhases.push_back(0);
    grid[nucleus[0]][nucleus[1]].grainPhases[numGrains - 1] = 1.0;
    orientations[numGrains] = {dist(gen),dist(gen),dist(gen)};
}
       
      