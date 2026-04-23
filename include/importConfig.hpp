#pragma once
#include "json.hpp" 
using json = nlohmann::json;
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>


struct config {

    double particleRadius; //m
    double particleLiquidIntEnergy; //J/m^2
    double particleSolidIntEnergy; //J/m^2

    double meltTemp; //K
    double dx; //um
    double dt; //s
    double kSolid; //W/m*k
    double kLiquid; //W/m*k
    int success;
    double startTemp; //K
    double cellWidth;//um
    double cellHeight;//um
    double cellArea;
    double particleVolFraction;
    double liqSolIntWidth; //um
    double grainIntWidth; //um
    int timeSteps;
    double underCoolReq;
    int totalSteps;
    double liqSolIntE;
    double grainIntE;
    std::vector<int> steps;
    double phaseCoefficient;
    double barrierHeightPhase;
    double barrierHeightGrain;
    double phasePreCo;
    double grainPreCo;
    double particleSlowingCoefficient;
    double phaseGradCo;
    double basePlateTemp;
    double heatCapacity;
    double density;
    double grainGradCo;
    double homoNucCoeff;
    double minTemp;
    double tGrad;
    double coolingRate;
    double drivingForceSlopek;
    double drivingForceIntercept;
    double molarMass;
    double molarVolume;
    double diffusionActivationEnergy;
    double particleDiameter; //m
    double hetNucUnderCooling;
    double particleDensity; //kg/m^3
    double thickness2D; //um

    bool radialCooling;                 // Enable radial cooling with colder outside
    bool hemisphericalCooling;          // Enable hemispherical temperature distribution hot at the top center
    
    // Grid dimensions
    int gridRows;     // Number of rows in computational grid
    int gridCols;     // Number of columns in computational grid
    
    // Configuration flags
    bool enableVisualization;  // Enable/disable live visualization
    bool enableProfiling;      // Enable/disable performance profiling

};


config inputConfig(std::string configSource);