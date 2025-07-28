#include <iostream>
#include "modelStartup.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <cmath>
config initModel() {
    config config{};
    std::string configPath = "../config/modelConfig.json";
    std::ifstream configIn(configPath);
    if (!configIn.is_open()) {
        std::cout << "Model Config File Not Found. Generating default config at: " << configPath << std::endl;
        std::ofstream configOut(configPath);
        if (!configOut) {
            std::cerr << "Failed to creat config file!" << std::endl;
            config.success =  0;
        }
        configOut << "{\n";
        configOut << "  \"dx_um\": 0.5,\n";
        configOut << "  \"dt_s\": 0.001,\n";
        configOut << "  \"startTempK\": 3000,\n";
        configOut << "  \"solidConductivityWMK\": 138,\n";
        configOut << "  \"cellHeightum\": 100,\n";
        configOut << "  \"cellWidthum\": 100\n";
        configOut << "}\n";
        configOut.close();
        std::cout << "Default Config Created. Please Edit and Rerun Program";
        config.success = 0;


    }
    std::string line;
    while (std::getline(configIn, line)) {
        std::istringstream iss(line);
        std::string key, eq;
        double value;

        if (!(iss >> key >> eq >> value)) continue;

        if (key == "dx_um") config.dx = value;
        if (key == "dt_s") config.dt = value;
        if (key == "startTempK") config.startTemp = value;
        if (key == "soludConductivityWMK") config.kSolid = value;
        if (key == "cellWidthum") config.cellWidth = value;
        if (key == "cellHeightum") config.cellHeight = value;
     
        std::cout << "Config Line: "<< line << std::endl;


    }

    configIn.close();





    fields starterField{};
    starterField.steps = {
    static_cast<int>(std::ceil(config.cellWidth / config.dx)),
    static_cast<int>(std::ceil(config.cellHeight / config.dx))
    };
    int totalSteps = starterField.steps[0] * starterField.steps[1];

    std::vector<double> tempTemp(totalSteps,config.startTemp);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0,25);

    for (double& T : tempTemp) {
        T += dist(gen);
    };
    starterField.temp = tempTemp;
    starterField.phase = std::vector<double>(totalSteps,0);
    starterField.orientation = {
        std::vector<double>(totalSteps,0), //Theta1
        std::vector<double>(totalSteps,0), //Phi
        std::vector<double>(totalSteps,0) //Theta2
    };






    config.success = 1;
    return config;
}