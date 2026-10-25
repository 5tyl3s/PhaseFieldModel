#include "json.hpp"
using json = nlohmann::json;

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>
#include "importConfig.hpp"

namespace {

double nucleationBarrier(
    double radius,
    double substrateRadius,
    double drivingForce,
    double gammaSolidSubstrate,
    double gammaLiquidSubstrate,
    double gammaLiquidSolid,
    double contactAngle) {
    const double substrateCapAngle = std::atan2(
        radius * std::sin(contactAngle),
        substrateRadius - radius * std::cos(contactAngle));
    const double nucleusCapAngle = std::atan2(
        substrateRadius * std::sin(contactAngle),
        substrateRadius * std::cos(contactAngle) - radius);

    const double substrateArea = 2.0 * M_PI * substrateRadius * substrateRadius *
        (1.0 - std::cos(substrateCapAngle));
    const double nucleusArea = 2.0 * M_PI * radius * radius *
        (1.0 - std::cos(nucleusCapAngle));
    const double capCosine = std::cos(nucleusCapAngle);
    const double capOffset = 1.0 - capCosine;
    const double nucleusVolume = (M_PI / 3.0) * radius * radius * radius *
        capOffset * capOffset * (2.0 + capCosine);

    return -nucleusVolume * drivingForce +
        substrateArea * (gammaSolidSubstrate - gammaLiquidSubstrate) +
        nucleusArea * gammaLiquidSolid;
}

bool peakNucleationBarrier(
    double substrateRadius,
    double drivingForce,
    double gammaSolidSubstrate,
    double gammaLiquidSubstrate,
    double gammaLiquidSolid,
    double contactAngle,
    double& peakBarrier) {
    const double classicalRadius = 2.0 * gammaLiquidSolid / drivingForce;
    const double minimumRadius = 1e-6 * std::min(substrateRadius, classicalRadius);
    const double maximumRadius = 30.0 * std::max(substrateRadius, classicalRadius);
    constexpr int radiusSamples = 1200;

    if (!std::isfinite(minimumRadius) || !std::isfinite(maximumRadius) ||
        minimumRadius <= 0.0 || maximumRadius <= minimumRadius) {
        return false;
    }

    std::vector<double> radii(radiusSamples);
    std::vector<double> barriers(radiusSamples);
    const double logMinimum = std::log(minimumRadius);
    const double logMaximum = std::log(maximumRadius);
    int peakIndex = 0;

    for (int i = 0; i < radiusSamples; ++i) {
        const double fraction = static_cast<double>(i) / (radiusSamples - 1);
        radii[i] = std::exp(logMinimum + fraction * (logMaximum - logMinimum));
        barriers[i] = nucleationBarrier(
            radii[i], substrateRadius, drivingForce, gammaSolidSubstrate,
            gammaLiquidSubstrate, gammaLiquidSolid, contactAngle);
        if (!std::isfinite(barriers[i])) {
            return false;
        }
        if (barriers[i] > barriers[peakIndex]) {
            peakIndex = i;
        }
    }

    if (peakIndex == 0 || peakIndex == radiusSamples - 1) {
        return false;
    }

    double lower = radii[peakIndex - 1];
    double upper = radii[peakIndex + 1];
    constexpr double goldenRatio = 0.6180339887498948482;
    double x1 = upper - goldenRatio * (upper - lower);
    double x2 = lower + goldenRatio * (upper - lower);
    double barrier1 = nucleationBarrier(
        x1, substrateRadius, drivingForce, gammaSolidSubstrate,
        gammaLiquidSubstrate, gammaLiquidSolid, contactAngle);
    double barrier2 = nucleationBarrier(
        x2, substrateRadius, drivingForce, gammaSolidSubstrate,
        gammaLiquidSubstrate, gammaLiquidSolid, contactAngle);

    for (int iteration = 0; iteration < 80; ++iteration) {
        if (barrier1 < barrier2) {
            lower = x1;
            x1 = x2;
            barrier1 = barrier2;
            x2 = lower + goldenRatio * (upper - lower);
            barrier2 = nucleationBarrier(
                x2, substrateRadius, drivingForce, gammaSolidSubstrate,
                gammaLiquidSubstrate, gammaLiquidSolid, contactAngle);
        } else {
            upper = x2;
            x2 = x1;
            barrier2 = barrier1;
            x1 = upper - goldenRatio * (upper - lower);
            barrier1 = nucleationBarrier(
                x1, substrateRadius, drivingForce, gammaSolidSubstrate,
                gammaLiquidSubstrate, gammaLiquidSolid, contactAngle);
        }
    }

    peakBarrier = std::max(barrier1, barrier2);
    return std::isfinite(peakBarrier);
}

bool heterogeneousNucleationThreshold(
    const config& modelConfig,
    double& thresholdTemperature,
    bool& completeWetting) {
    const double gammaSolidSubstrate = modelConfig.particleSolidIntEnergy;
    const double gammaLiquidSubstrate = modelConfig.particleLiquidIntEnergy;
    const double gammaLiquidSolid = modelConfig.liqSolIntE;
    const double substrateRadius = modelConfig.particleRadius;
    const double minimumUndercooling =
        std::max(1e-3, std::numeric_limits<double>::epsilon() * modelConfig.meltTemp);

    if (!std::isfinite(substrateRadius) || substrateRadius <= 0.0 ||
        !std::isfinite(gammaSolidSubstrate) || gammaSolidSubstrate <= 0.0 ||
        !std::isfinite(gammaLiquidSubstrate) || gammaLiquidSubstrate <= 0.0 ||
        !std::isfinite(gammaLiquidSolid) || gammaLiquidSolid <= 0.0 ||
        !std::isfinite(modelConfig.meltTemp) || modelConfig.meltTemp <= minimumUndercooling ||
        !std::isfinite(modelConfig.molarVolume) || modelConfig.molarVolume <= 0.0 ||
        !std::isfinite(modelConfig.drivingForceSlopek) ||
        !std::isfinite(modelConfig.drivingForceIntercept) ||
        modelConfig.drivingForceSlopek == 0.0 ||
        !std::isfinite(modelConfig.underCoolReq) || modelConfig.underCoolReq <= 0.0) {
        std::cerr << "[config] Cannot calculate heterogeneous nucleation threshold: "
                     "invalid thermodynamic or particle parameters." << std::endl;
        return false;
    }

    const double rawCosine =
        (gammaLiquidSubstrate - gammaSolidSubstrate) / gammaLiquidSolid;
    if (!std::isfinite(rawCosine)) {
        std::cerr << "[config] Cannot calculate heterogeneous nucleation threshold: "
                     "invalid Young contact-angle cosine." << std::endl;
        return false;
    }

    completeWetting = std::abs(rawCosine) >= 1.0;
    if (completeWetting) {
        std::cerr << "[config] Young contact-angle cosine is outside the partial-wetting "
                     "range; using the complete-wetting limit." << std::endl;
        thresholdTemperature = modelConfig.meltTemp;
        return true;
    }

    const double contactAngle = std::acos(rawCosine);
    const auto drivingForceAt = [&modelConfig](double temperature) {
        return -(modelConfig.drivingForceSlopek * temperature +
                 modelConfig.drivingForceIntercept) / modelConfig.molarVolume;
    };

    const double referenceTemperature = modelConfig.meltTemp - modelConfig.underCoolReq;
    const double referenceDrivingForce = drivingForceAt(referenceTemperature);
    if (!std::isfinite(referenceDrivingForce) || referenceDrivingForce <= 0.0) {
        std::cerr << "[config] Cannot calculate heterogeneous nucleation threshold: "
                     "reference undercooling does not produce a positive driving force."
                  << std::endl;
        return false;
    }

    const double referenceHomogeneousBarrier =
        16.0 * M_PI * std::pow(gammaLiquidSolid, 3.0) /
        (3.0 * referenceDrivingForce * referenceDrivingForce);
    const double targetBarrier = 1e-3 * referenceHomogeneousBarrier;
    const double minimumTemperature = modelConfig.minTemp;
    const double maximumTemperature = modelConfig.meltTemp - minimumUndercooling;

    if (!std::isfinite(minimumTemperature) ||
        minimumTemperature >= maximumTemperature) {
        std::cerr << "[config] Cannot calculate heterogeneous nucleation threshold: "
                     "temperature search interval is empty." << std::endl;
        return false;
    }

    double lowBarrier = 0.0;
    double highBarrier = 0.0;
    if (!peakNucleationBarrier(
            substrateRadius, drivingForceAt(minimumTemperature),
            gammaSolidSubstrate, gammaLiquidSubstrate, gammaLiquidSolid,
            contactAngle, lowBarrier) ||
        !peakNucleationBarrier(
            substrateRadius, drivingForceAt(maximumTemperature),
            gammaSolidSubstrate, gammaLiquidSubstrate, gammaLiquidSolid,
            contactAngle, highBarrier)) {
        std::cerr << "[config] Cannot calculate heterogeneous nucleation threshold: "
                     "a resolved interior barrier peak was not found." << std::endl;
        return false;
    }

    double lowerTemperature = minimumTemperature;
    double upperTemperature = maximumTemperature;
    double lowerDifference = lowBarrier - targetBarrier;
    double upperDifference = highBarrier - targetBarrier;
    if (lowerDifference == 0.0) {
        thresholdTemperature = lowerTemperature;
        return true;
    }
    if (upperDifference == 0.0) {
        thresholdTemperature = upperTemperature;
        return true;
    }
    if ((lowerDifference < 0.0) == (upperDifference < 0.0)) {
        std::cerr << "[config] Cannot calculate heterogeneous nucleation threshold: "
                     "target barrier is not crossed in the configured temperature range."
                  << std::endl;
        return false;
    }

    for (int iteration = 0; iteration < 80; ++iteration) {
        const double midpoint = 0.5 * (lowerTemperature + upperTemperature);
        double midpointBarrier = 0.0;
        if (!peakNucleationBarrier(
                substrateRadius, drivingForceAt(midpoint),
                gammaSolidSubstrate, gammaLiquidSubstrate, gammaLiquidSolid,
                contactAngle, midpointBarrier)) {
            std::cerr << "[config] Cannot calculate heterogeneous nucleation threshold: "
                         "a resolved interior barrier peak was not found." << std::endl;
            return false;
        }

        const double midpointDifference = midpointBarrier - targetBarrier;
        if ((midpointDifference < 0.0) == (lowerDifference < 0.0)) {
            lowerTemperature = midpoint;
            lowerDifference = midpointDifference;
        } else {
            upperTemperature = midpoint;
            upperDifference = midpointDifference;
        }
    }

    thresholdTemperature = 0.5 * (lowerTemperature + upperTemperature);
    return std::isfinite(thresholdTemperature);
}

// Tabulate the peak heterogeneous barrier versus temperature so the
// simulation can evaluate the nucleation rate cheaply every time step.
void buildHetBarrierTable(config& modelConfig) {
    constexpr int tableSize = 600;
    const double inf = std::numeric_limits<double>::infinity();
    const double minimumUndercooling =
        std::max(1e-3, std::numeric_limits<double>::epsilon() * modelConfig.meltTemp);
    const double tMin = modelConfig.minTemp;
    const double tMax = modelConfig.meltTemp - minimumUndercooling;

    modelConfig.hetBarrierTable.clear();
    if (!(tMax > tMin) || modelConfig.molarVolume <= 0.0 ||
        modelConfig.liqSolIntE <= 0.0 || modelConfig.particleRadius <= 0.0) {
        return;
    }

    const double rawCosine =
        (modelConfig.particleLiquidIntEnergy - modelConfig.particleSolidIntEnergy) /
        modelConfig.liqSolIntE;
    const bool completeWetting = !std::isfinite(rawCosine) || std::abs(rawCosine) >= 1.0;
    const double contactAngle = completeWetting ? 0.0 : std::acos(rawCosine);

    modelConfig.hetBarrierTempMin = tMin;
    modelConfig.hetBarrierTempStep = (tMax - tMin) / (tableSize - 1);
    modelConfig.hetBarrierTable.assign(tableSize, inf);

    for (int i = 0; i < tableSize; ++i) {
        if (completeWetting) {
            modelConfig.hetBarrierTable[i] = 0.0;
            continue;
        }
        const double temperature = tMin + i * modelConfig.hetBarrierTempStep;
        const double drivingForce = -(modelConfig.drivingForceSlopek * temperature +
                                      modelConfig.drivingForceIntercept) / modelConfig.molarVolume;
        double barrier = 0.0;
        if (drivingForce > 0.0 && peakNucleationBarrier(
                modelConfig.particleRadius, drivingForce,
                modelConfig.particleSolidIntEnergy, modelConfig.particleLiquidIntEnergy,
                modelConfig.liqSolIntE, contactAngle, barrier) && barrier >= 0.0) {
            modelConfig.hetBarrierTable[i] = barrier;
        }
    }
}

}

double hetNucPeakBarrier(double temp, const config& modelConf) {
    const auto& table = modelConf.hetBarrierTable;
    const double inf = std::numeric_limits<double>::infinity();
    if (table.size() < 2 || modelConf.hetBarrierTempStep <= 0.0) return inf;

    const double pos = (temp - modelConf.hetBarrierTempMin) / modelConf.hetBarrierTempStep;
    if (pos < 0.0 || pos > static_cast<double>(table.size() - 1)) return inf;

    size_t i0 = static_cast<size_t>(pos);
    if (i0 >= table.size() - 1) i0 = table.size() - 2;
    const double b0 = table[i0];
    const double b1 = table[i0 + 1];
    if (!std::isfinite(b0) || !std::isfinite(b1)) return inf;
    const double f = pos - static_cast<double>(i0);
    return b0 + f * (b1 - b0);
}

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
        configOut << "    \"CoolingRatempers\": \"Cooling front speed (m/s).\",\n";
        configOut << "    \"tempGradKperm\": \"Vertical temperature gradient (K per meter).\",\n";
        configOut << "    \"radialCooling\": \"Enable radial cooling with colder outside.\",\n";
        configOut << "    \"hemisphericalCooling\": \"Enable hemispherical temperature distribution hot at the top center.\",\n";
        configOut << "    \"periodicLeft\": \"Enable periodic boundary conditions on left edge.\",\n";
        configOut << "    \"periodicRight\": \"Enable periodic boundary conditions on right edge.\",\n";
        configOut << "    \"periodicTop\": \"Enable periodic boundary conditions on top edge.\",\n";
        configOut << "    \"periodicBottom\": \"Enable periodic boundary conditions on bottom edge.\",\n";
        configOut << "    \"simulationTimeLimit_s\": \"Maximum wall-clock runtime in seconds (0 = disabled).\",\n";
        configOut << "    \"MinimumTempK\": \"Minimum allowed temperature (K).\",\n";
        configOut << "    \"molarMass\": \"Molar mass of material (g/mol).\",\n";
        configOut << "    \"drivingForceSlopek\": \"Slope used in driving force calculation.\",\n";
        configOut << "    \"drivingForceIntercept\": \"Intercept used in driving force calculation.\",\n";
        configOut << "    \"particleSolidIntEnergy\": \"Particle interaction energy in solid (model units).\",\n";
        configOut << "    \"particleLiqIntEnergy\": \"Particle interaction energy in liquid (model units).\",\n";
        configOut << "    \"PhaseMobility\": \"Phase mobility coefficient used in the phase evolution update.\",\n";
        configOut << "    \"GrainMobility\": \"Grain mobility coefficient used in the grain evolution update.\",\n";
        configOut << "    \"ParticleMobilityLiquid\": \"Particle mobility in liquid phase.\",\n";
        configOut << "    \"ParticleMobilitySolid\": \"Particle mobility in solid phase.\"\n";
        configOut << "  },\n";
        configOut << "  \"dx_um\": 1,\n";
        configOut << "  \"thickness2D_um\": 1,\n";
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
        configOut << "  \"CoolingRatempers\": 0.001,\n";
        configOut << "  \"MinimumTempK\": 100,\n";
        configOut << "  \"tempGradKperm\": 20000,\n";
        configOut << "  \"radialCooling\": false,\n";
        configOut << "  \"hemisphericalCooling\": false,\n";
        configOut << "  \"periodicLeft\": true,\n";
        configOut << "  \"periodicRight\": true,\n";
        configOut << "  \"periodicTop\": false,\n";
        configOut << "  \"periodicBottom\": false,\n";
        configOut << "  \"simulationTimeLimit_s\": 0.0,\n";
        configOut << "  \"molarMass\": 95.95,\n";
        configOut << "  \"drivingForceSlopek\": 13.8,\n";
        configOut << "  \"drivingForceIntercept\": 39842,\n";
        configOut << "  \"particleSolidIntEnergy\": 2.0,\n";
        configOut << "  \"particleLiqIntEnergy\": 1.0,\n";
        configOut << "  \"PhaseMobility\": 4e-14,\n";
        configOut << "  \"GrainMobility\": 2e-14,\n";
        configOut << "  \"ParticleMobilityLiquid\": 2e-4,\n";
        configOut << "  \"ParticleMobilitySolid\": 1e-18,\n";
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
    newConfig.thickness2D = readOr<double>(j, "thickness2D_um", 1.0) * 1e-6;
    newConfig.grainIntE = readOr<double>(j, "GrainIntEnergyJm2", 3.0);
    newConfig.grainIntWidth = readOr<double>(j, "GrainIntWidthum", 1.0) * 1e-6;

    newConfig.barrierHeightPhase = readOr<double>(j, "PhaseBarrierHeightCoefficient", 0.25);
    newConfig.barrierHeightGrain = readOr<double>(j, "GrainBarrierHeightCoefficient", 0.125);

    newConfig.phaseCoefficient = 0.75  * newConfig.liqSolIntE / (newConfig.liqSolIntWidth*newConfig.barrierHeightPhase);
    newConfig.phaseMobility = readOr<double>(j, "PhaseMobility", 4e-14);
    newConfig.grainMobility = readOr<double>(j, "GrainMobility", 2e-14);
    newConfig.particleMobilityLiquid = readOr<double>(j, "ParticleMobilityLiquid", 2e-4);
    newConfig.particleMobilitySolid = readOr<double>(j, "ParticleMobilitySolid", 1e-18);
    
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

    newConfig.drivingForceSlopek = readOr<double>(j, "drivingForceSlopek", 13.8)*1000;
    newConfig.drivingForceIntercept = readOr<double>(j, "drivingForceIntercept", 39842)*1000;
    newConfig.diffusionActivationEnergy = readOr<double>(j, "diffusionActivationEnergy", 1.5);

    newConfig.molarMass = readOr<double>(j, "molarMass", 95.95);
    newConfig.molarVolume = newConfig.molarMass / (newConfig.density * 1000);
    std::cout << "Molar Volume: " << newConfig.molarVolume << " m^3/mol" << std::endl;

    newConfig.particleSolidIntEnergy = readOr<double>(j, "particleSolidIntEnergy", 2.0);
    newConfig.particleLiquidIntEnergy = readOr<double>(j, "particleLiqIntEnergy", 1.0);
    newConfig.particleDiameter = readOr<double>(j, "particleDiameter_um", 1.0) * 1e-6;

    newConfig.radialCooling = readOr<bool>(j, "radialCooling", false);
    newConfig.hemisphericalCooling = readOr<bool>(j, "hemisphericalCooling", false);

    // Periodic boundary conditions
    newConfig.periodicLeft = readOr<bool>(j, "periodicLeft", true);
    newConfig.periodicRight = readOr<bool>(j, "periodicRight", true);
    newConfig.periodicTop = readOr<bool>(j, "periodicTop", false);
    newConfig.periodicBottom = readOr<bool>(j, "periodicBottom", false);

    std::cout << "Particle Radius set to: " << newConfig.particleRadius << " m" << std::endl;
    // Discretization steps
    int height = static_cast<int>(std::ceil(newConfig.cellHeight / newConfig.dx));
    int width  = static_cast<int>(std::ceil(newConfig.cellWidth / newConfig.dx));
    newConfig.steps = {height, width};
    newConfig.totalSteps = height * width;
    double energyTempConv = (newConfig.heatCapacity)/(newConfig.density)*newConfig.dx*newConfig.dx; // J/K
    // Match the MATLAB barrier criterion: threshold where the peak heterogeneous
    // barrier reaches 0.1% of the homogeneous barrier at the reference undercooling.
    bool completeWetting = false;
    if (!heterogeneousNucleationThreshold(
            newConfig, newConfig.hetNucUnderCooling, completeWetting)) {
        return newConfig;
    }
    std::cout << "Heterogeneous nucleation threshold: "
              << newConfig.hetNucUnderCooling << " K"
              << (completeWetting ? " (complete-wetting limit)" : "") << std::endl;
    std::cout << "Heterogeneous Nucleation Undercooling Threshold: " << newConfig.hetNucUnderCooling << " K" << std::endl;
    // Sites per particle = 4*pi*Rs^2 / a^2, with a the mean atomic spacing from the molar volume.
    // n0 = N_p * N_site/particle is then particleCount * sites per cell in the nucleation check.
    {
        const double avogadro = 6.02214076e23;
        const double atomicSpacing = std::cbrt(newConfig.molarVolume / avogadro);
        const double defaultSites = 4.0 * M_PI * newConfig.particleRadius * newConfig.particleRadius /
                                    (atomicSpacing * atomicSpacing);
        newConfig.hetNucSitesPerParticle = readOr<double>(j, "hetNucSitesPerParticle", defaultSites);
    }
    buildHetBarrierTable(newConfig);
    // Configuration flags
    newConfig.enableVisualization = readOr<bool>(j, "enableVisualization", true);
    newConfig.enableProfiling = readOr<bool>(j, "enableProfiling", true);
    newConfig.simulationTimeLimit = readOr<double>(j, "simulationTimeLimit_s", 0.0);
    
    // Grid dimensions
    newConfig.gridRows = readOr<int>(j, "gridRows", 500);
    newConfig.gridCols = readOr<int>(j, "gridCols", 500);
    
    newConfig.particleDensity = 3150;
    newConfig.success = 1; // successful load
    return newConfig;
}
