#include "json.hpp" 
using json = nlohmann::json;
#include <iostream>

#include "importConfig.hpp"
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

config inputConfig(std::string configSource) {
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
        configOut << "  \"HeatCapacityJKgK\": 251,\n"; // Mo
        configOut << "  \"Densitykgm3\": 10280,\n";    // Mo
        configOut << "  \"dt_s\": 0.001,\n";
        configOut << "  \"startTempK\": 3000,\n";
        configOut << "  \"solidConductivityWMK\": 138,\n"; // Mo
        configOut << "  \"cellHeightum\": 100,\n";
        configOut << "  \"cellWidthum\": 100,\n";
        configOut << "  \"particleVolumeFraction\": 0.01,\n";
        configOut << "  \"LiqSolIntEnergyJm2\": 0.3,\n";
        configOut << "  \"LiqSolIntWidthum\": 0.3,\n";
        configOut << "  \"GrainIntEnergyJm2\": 0.3,\n";
        configOut << "  \"GrainIntWidthum\": 0.3,\n";
        configOut << "  \"PhaseBarrierHeightCoefficient\": 0.25,\n";
        configOut << "  \"BasePlateTempK\": 500.00,\n";
        configOut << "  \"HomogeneousNucleationCoefficient\": 100,\n";
        configOut << "  \"TimeSteps\": 10000,\n";
        configOut << "  \"UndercoolingRequirement\": 0.95,\n";
        configOut << "  \"CoolingRatempers\": 100,\n";
        configOut << "  \"MinimumTempK\": 300,\n";
        configOut << "  \"tempGradKperm\": 1.0,\n";
        configOut << "  \"molarMass\": 95.95,\n";
        configOut << "  \"drivingForceSlopek\": 0.0138,\n";
        configOut << "  \"drivingForceIntercept\": 39.842,\n";
        configOut << "  \"diffusionActivationEnergy\": 1.5,\n";
        configOut << "  \"ParticleSolidIntEnergy\": 6.0,\n";
        configOut << "  \"ParticleLiqIntEnergy\": 3.0,\n";
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
    newConfig.liqSolIntWidth = newConfig.liqSolIntWidth*1e-6;
    newConfig.phaseCoefficient = 0.75*newConfig.liqSolIntWidth * newConfig.liqSolIntE;
    newConfig.barrierHeightGrain = j["GrainBarrierHeightCoefficient"];
    newConfig.barrierHeightGrain = newConfig.barrierHeightGrain*1e-6;
    newConfig.barrierHeightPhase = j["PhaseBarrierHeightCoefficient"];
    newConfig.barrierHeightPhase = newConfig.barrierHeightPhase*1e-6;
    newConfig.basePlateTemp = j["BasePlateTempK"];
    newConfig.meltTemp = j["meltTemp"];
    newConfig.heatCapacity = j["HeatCapacityJKgK"];
    newConfig.density = j["Densitykgm3"];
    newConfig.kLiquid = newConfig.kSolid * 0.45; 
    newConfig.grainIntWidth = j["GrainIntWidthum"];
    newConfig.grainIntWidth = newConfig.grainIntWidth*1e-6;
    newConfig.grainIntE = j["GrainIntEnergyJm2"];
    newConfig.grainIntE = newConfig.grainIntE;
    newConfig.phasePreCo = 0.75*newConfig.liqSolIntE/(newConfig.barrierHeightPhase*newConfig.liqSolIntWidth);
    newConfig.grainPreCo = 0.75*newConfig.grainIntE/(newConfig.barrierHeightGrain*newConfig.grainIntWidth);
    newConfig.phaseGradCo = 0.75*newConfig.liqSolIntE*newConfig.liqSolIntWidth;
    newConfig.particleSlowingCoefficient = 0;
    newConfig.success = 1;
    newConfig.timeSteps = j["TimeSteps"];
    newConfig.grainGradCo = 0.5*newConfig.grainIntWidth;
    newConfig.homoNucCoeff =j["HomogeneousNucleationCoefficient"];
    newConfig.underCoolReq = j["UndercoolingRequirement"];
    newConfig.coolingRate = j["CoolingRatempers"];
    std::cout << newConfig.coolingRate << std::endl;
    newConfig.tGrad = j["tempGradKperm"];
    newConfig.coolingRate =j["CoolingRatempers"];
    newConfig.minTemp = j["MinimumTempK"];
    newConfig.drivingForceSlopek = j["drivingForceSlopek"];
    newConfig.drivingForceIntercept = j["drivingForceIntercept"];
    newConfig.diffusionActivationEnergy = j["diffusionActivationEnergy"];
    newConfig.molarMass = j["molarMass"];
    newConfig.particleSolidIntEnergy = j["particleSolidIntEnergy"];
    newConfig.particleLiquidIntEnergy = j["particleLiqIntEnergy"];

    newConfig.molarVolume = (newConfig.molarMass) / (newConfig.density*1000);

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