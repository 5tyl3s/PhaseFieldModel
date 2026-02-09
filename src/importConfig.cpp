#include "json.hpp"
using json = nlohmann::json;

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "importConfig.hpp"

// Helper: safely read JSON value or fallback to default
template<typename T>
T readOr(const json& j, const std::string& key, T defaultValue) {
    if (j.contains(key) && !j[key].is_null()) {
        return j[key].get<T>();
    } else {
        std::cerr << "[config] Key '" << key << "' missing or null. Using default: " << defaultValue << std::endl;
        return defaultValue;
    }
}

config inputConfig(std::string configSource) {
    config newConfig;
    newConfig.success = 0; // pessimistic until finished

    std::string configPath = configSource.empty() ? "../config/modelConfig.json" : configSource;
    std::ifstream configIn(configPath);

    // If file does not exist, create a default JSON
    if (!configIn.is_open()) {
        std::cout << "[config] Config file not found. Generating default at " << configPath << std::endl;
        std::ofstream configOut(configPath);
        if (!configOut) {
            std::cerr << "[config] Failed to create config file!" << std::endl;
            return newConfig;
        }

        // Write default values matching your provided constants (reorganized, with descriptions stored separately in config file)
        configOut << "{\n";
        configOut << "  \"__comments\": {\n";
        configOut << "    \"dx_um\": \"Grid spacing in micrometers (um). Converted to meters in code: dx = dx_um * 1e-6\",\n";
        configOut << "    \"dt_s\": \"Time step in seconds.\",\n";
        configOut << "    \"startTempK\": \"Initial temperature in Kelvin.\",\n";
        configOut << "    \"meltTemp\": \"Material melting temperature in Kelvin.\",\n";
        configOut << "    \"particleVolumeFraction\": \"Initial particle volume fraction (0..1).\",\n";
        configOut << "    \"particleDiameter_um\": \"Particle diameter in micrometers (um).\",\n";
        configOut << "    \"LiqSolIntEnergyJm2\": \"Liquid–solid interfacial energy (J/m^2).\",\n";
        configOut << "    \"LiqSolIntWidthum\": \"Liquid–solid interfacial width (um).\",\n";
        configOut << "    \"GrainIntEnergyJm2\": \"Grain boundary energy (J/m^2).\",\n";
        configOut << "    \"GrainIntWidthum\": \"Grain boundary width (um).\",\n";
        configOut << "    \"PhaseBarrierHeightCoefficient\": \"Phase-field barrier height coefficient (dimensionless).\",\n";
        configOut << "    \"GrainBarrierHeightCoefficient\": \"Grain-field barrier height coefficient (dimensionless).\",\n";
        configOut << "    \"BasePlateTempK\": \"Base plate temperature (K) used for boundary/initial conditions.\",\n";
        configOut << "    \"diffusionActivationEnergy\": \"Activation energy for diffusion (units as used in nucleation formula).\",\n";
        configOut << "    \"TimeSteps\": \"Number of timesteps to run the simulation.\",\n";
        configOut << "    \"CoolingRatempers\": \"Cooling rate (K per second).\",\n";
        configOut << "    \"tempGradKperm\": \"Vertical temperature gradient (K per meter).\",\n";
        configOut << "    \"MinimumTempK\": \"Minimum allowed temperature (K).\",\n";
        configOut << "    \"molarMass\": \"Molar mass of material (g/mol).\",\n";
        configOut << "    \"drivingForceSlopek\": \"Slope used in driving force calculation.\",\n";
        configOut << "    \"drivingForceIntercept\": \"Intercept used in driving force calculation.\",\n";
        configOut << "    \"particleSolidIntEnergy\": \"Particle interaction energy in solid (model units).\",\n";
        configOut << "    \"particleLiqIntEnergy\": \"Particle interaction energy in liquid (model units).\"\n";
        configOut << "  },\n";
        configOut << "  \"dx_um\": 1,\n";
        configOut << "  \"dt_s\": 1e-12,\n";
        configOut << "  \"startTempK\": 2896,\n";
        configOut << "  \"meltTemp\": 2896,\n";
        configOut << "  \"HeatCapacityJKgK\": 251,\n";
        configOut << "  \"Densitykgm3\": 10280,\n";
        configOut << "  \"cellHeightum\": 50,\n";
        configOut << "  \"cellWidthum\": 50,\n";
        configOut << "  \"particleVolumeFraction\": 0.0,\n";
        configOut << "  \"particleDiameter_um\": 1.0,\n";
        configOut << "  \"LiqSolIntEnergyJm2\": 3,\n";
        configOut << "  \"LiqSolIntWidthum\": 1,\n";
        configOut << "  \"GrainIntEnergyJm2\": 3,\n";
        configOut << "  \"GrainIntWidthum\": 1,\n";
        configOut << "  \"PhaseBarrierHeightCoefficient\": 0.25,\n";
        configOut << "  \"BasePlateTempK\": 500.0,\n";
        configOut << "  \"diffusionActivationEnergy\": 1.5,\n";
        configOut << "  \"TimeSteps\": 1000,\n";
        configOut << "  \"CoolingRatempers\": 3e11,\n";
        configOut << "  \"MinimumTempK\": 100,\n";
        configOut << "  \"tempGradKperm\": 20000,\n";
        configOut << "  \"molarMass\": 95.95,\n";
        configOut << "  \"drivingForceSlopek\": 13.8,\n";
        configOut << "  \"drivingForceIntercept\": 39842,\n";
        configOut << "  \"particleSolidIntEnergy\": 2.0,\n";
        configOut << "  \"particleLiqIntEnergy\": 1.0,\n";
        configOut << "  \"GrainBarrierHeightCoefficient\": 0.125,\n";
        configOut << "  \"enableVisualization\": true,\n";
        configOut << "  \"enableProfiling\": true\n";
        configOut << "}\n";
        configOut.close();
        std::cout << "[config] Default config created. Please edit and re-run program.\n";
    }

    // Read JSON
    json j;
    configIn.open(configPath);
    if (configIn.is_open()) {
        configIn >> j;
    }

    // -----------------------------
    // Read values with defaults
    // -----------------------------
    newConfig.dx = readOr<double>(j, "dx_um", 1.0) * 1e-6;
    newConfig.dt = readOr<double>(j, "dt_s", 1e-3);
    newConfig.startTemp = readOr<double>(j, "startTempK", 2896);
    newConfig.kSolid = readOr<double>(j, "solidConductivityWMK", 138);
    newConfig.cellHeight = readOr<double>(j, "cellHeightum", 50) * 1e-6;
    newConfig.cellWidth = readOr<double>(j, "cellWidthum", 50) * 1e-6;
    newConfig.cellArea = newConfig.cellHeight * newConfig.cellWidth;
    newConfig.particleVolFraction = readOr<double>(j, "particleVolumeFraction", 0.0);
    // particle size (diameter) in micrometers -> store as meters
    newConfig.particleDiameter = readOr<double>(j, "particleDiameter_um", 1.0) * 1e-6;
    // also set particleRadius (legacy field) for older code
    newConfig.particleRadius = newConfig.particleDiameter * 0.5;

    newConfig.liqSolIntE = readOr<double>(j, "LiqSolIntEnergyJm2", 3.0);
    newConfig.liqSolIntWidth = readOr<double>(j, "LiqSolIntWidthum", 1.0) * 1e-6;

    newConfig.grainIntE = readOr<double>(j, "GrainIntEnergyJm2", 3.0);
    newConfig.grainIntWidth = readOr<double>(j, "GrainIntWidthum", 1.0) * 1e-6;

    newConfig.barrierHeightPhase = readOr<double>(j, "PhaseBarrierHeightCoefficient", 0.25);
    newConfig.barrierHeightGrain = readOr<double>(j, "GrainBarrierHeightCoefficient", 0.125);

    newConfig.phaseCoefficient = 0.75  * newConfig.liqSolIntE / (newConfig.liqSolIntWidth*newConfig.barrierHeightPhase);
    
    newConfig.phaseGradCo = 0.75 * newConfig.liqSolIntE * newConfig.liqSolIntWidth;
    newConfig.grainGradCo = 0.5 * newConfig.grainIntWidth;
    std::cout << "Grain Gradient Coefficient: " << newConfig.grainGradCo << std::endl;

    newConfig.phasePreCo = 0.75 * newConfig.liqSolIntE / (newConfig.barrierHeightPhase * newConfig.liqSolIntWidth);
    newConfig.grainPreCo = 0.75 * newConfig.grainIntE / (newConfig.barrierHeightGrain * newConfig.grainIntWidth);

    newConfig.basePlateTemp = readOr<double>(j, "BasePlateTempK", 500.0);
    newConfig.meltTemp = readOr<double>(j, "meltTemp", 2896);
    newConfig.heatCapacity = readOr<double>(j, "HeatCapacityJKgK", 251);
    newConfig.density = readOr<double>(j, "Densitykgm3", 10280);
    newConfig.kLiquid = newConfig.kSolid * 0.45;

    newConfig.particleSlowingCoefficient = 0;

    newConfig.timeSteps = readOr<int>(j, "TimeSteps", 1000);
    newConfig.homoNucCoeff = readOr<double>(j, "HomogeneousNucleationCoefficient", 1e20);
    newConfig.underCoolReq = readOr<double>(j, "UndercoolingRequirement", 1);
    newConfig.coolingRate = readOr<double>(j, "CoolingRatempers", 0.3);
    newConfig.tGrad = readOr<double>(j, "tempGradKperm", 20000);
    newConfig.minTemp = readOr<double>(j, "MinimumTempK", 100);

    newConfig.drivingForceSlopek = readOr<double>(j, "drivingForceSlopek", 13.8);
    newConfig.drivingForceIntercept = readOr<double>(j, "drivingForceIntercept", 39842);
    newConfig.diffusionActivationEnergy = readOr<double>(j, "diffusionActivationEnergy", 1.5);

    newConfig.molarMass = readOr<double>(j, "molarMass", 95.95);
    newConfig.molarVolume = newConfig.molarMass / (newConfig.density * 1000);
    std::cout << "Molar Volume: " << newConfig.molarVolume << " m^3/mol" << std::endl;

    newConfig.particleSolidIntEnergy = readOr<double>(j, "particleSolidIntEnergy", 2.0);
    newConfig.particleLiquidIntEnergy = readOr<double>(j, "particleLiqIntEnergy", 1.0);
    newConfig.particleDiameter = readOr<double>(j, "particleDiameter_um", 1.0) * 1e-6;

    std::cout << "Particle Radius set to: " << newConfig.particleRadius << " m" << std::endl;
    // Discretization steps
    int height = static_cast<int>(std::ceil(newConfig.cellHeight / newConfig.dx));
    int width  = static_cast<int>(std::ceil(newConfig.cellWidth / newConfig.dx));
    newConfig.steps = {height, width};
    newConfig.totalSteps = height * width;
    double energyTempConv = (newConfig.heatCapacity)/(newConfig.density)*newConfig.dx*newConfig.dx; // J/K
    // Match MATLAB: ((pi*r^2)*(Es-El) + (pi*(r+0.3e-9)^2)*liqSolIntE) / energyTempConv
    newConfig.hetNucUnderCooling = (
        (3.141592653589793 * pow(newConfig.particleRadius,2) * (newConfig.particleSolidIntEnergy - newConfig.particleLiquidIntEnergy))
        + (3.141592653589793 * pow((newConfig.particleRadius + 0.3e-9),2) * newConfig.liqSolIntE)
    ) / energyTempConv; // temperature change (K)
    std::cout << "Heterogeneous Nucleation Undercooling Threshold: " << newConfig.hetNucUnderCooling << " J/m^3" << std::endl;
    // Configuration flags
    newConfig.enableVisualization = readOr<bool>(j, "enableVisualization", true);
    newConfig.enableProfiling = readOr<bool>(j, "enableProfiling", true);
    newConfig.particleDensity = 3150;
    newConfig.success = 1; // successful load
    return newConfig;
}
