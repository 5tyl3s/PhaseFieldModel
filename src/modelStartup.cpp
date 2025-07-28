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
        configOut << "  \"dt_s\": 0.001,\n";
        configOut << "  \"startTempK\": 3000,\n";
        configOut << "  \"solidConductivityWMK\": 138,\n";
        configOut << "  \"cellHeightum\": 100,\n";
        configOut << "  \"cellWidthum\": 100,\n";
        configOut << "  \"particleVolumeFraction\": 0.01\n";
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

    for (int i = 0; i< heightSteps;i++) {
        for (int j = 0; j < widthSteps;j++) {
            node& gp = grid[i][j];
            gp.neighbors[0] = (j>0) ? &grid[i][j-1]:&grid[i][widthSteps-1]; //Left
            gp.neighbors[1] = (j<widthSteps-1) ? &grid[i][j+1]:&grid[i][0]; //Right
            gp.neighbors[2] = (i>0) ? &grid[i-1][j]: nullptr;
            gp.neighbors[3] = (i<heightSteps-1) ? &grid[i+1][j]: nullptr; //Down
        };
    };

};


void gridField::init(config modelConfig) {
    

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0,25);
    buildGrid(modelConfig.steps[0],modelConfig.steps[1]);

    for (int i = 0; i < modelConfig.steps[0]; i++) {
        for (int j = 0; j < modelConfig.steps[1]; j++) {
            grid[i][j].temp = modelConfig.startTemp +dist(gen);
            grid[i][j].phase = 0;
            grid[i][j].particleComp = modelConfig.particleVolFraction; 
        };
    }
}
      