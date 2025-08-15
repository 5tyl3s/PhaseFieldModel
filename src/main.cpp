#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <thread>

// These functions must be implemented to fill only the specified region of the output array
void calcGrainDiffEnergyRegion(gridField& model, config& configData,
    std::vector<std::vector<std::vector<double>>>& grainDiffEn,
    int i_start, int i_end, int j_start, int j_end);

void calcPhaseDiffEnergyRegion(gridField& model, config& configData,
    std::vector<std::vector<double>>& phaseDiffEn,
    int i_start, int i_end, int j_start, int j_end);

void calcTempDiffRegion(gridField& model, config& configData,
    std::vector<std::vector<double>>& tempGrad,
    int i_start, int i_end, int j_start, int j_end);

int main() {
    std::cout << "PhaseField Sim Starting" << std::endl;
    config configData = inputConfig();
    if (configData.success == 0) {
        std::cerr << "Config Failed";
        return 1;
    }

    gridField model;
    model.init(configData);
    std::vector<std::vector<std::vector<double>>> grainDiffEn(configData.steps[0],std::vector<std::vector<double>>(configData.steps[1],std::vector<double>(model.numGrains)));
    std::vector<std::vector<double>> phaseDiffEn(configData.steps[0],std::vector<double>(configData.steps[1]));
    std::vector<std::vector<double>> tempGrad(configData.steps[0],std::vector<double>(configData.steps[1]));

    int i_mid = configData.steps[0] / 2;
    int j_mid = configData.steps[1] / 2;

    for (int t = 0; t < 30000; t++) {
        std::cout << t << std::endl;

        // Launch threads to calculate each quadrant of each field
        std::thread th1([&](){
            calcGrainDiffEnergyRegion(model, configData, grainDiffEn, 0, i_mid, 0, j_mid);
            calcPhaseDiffEnergyRegion(model, configData, phaseDiffEn, 0, i_mid, 0, j_mid);
            calcTempDiffRegion(model, configData, tempGrad, 0, i_mid, 0, j_mid);
        });
        std::thread th2([&](){
            calcGrainDiffEnergyRegion(model, configData, grainDiffEn, 0, i_mid, j_mid, configData.steps[1]);
            calcPhaseDiffEnergyRegion(model, configData, phaseDiffEn, 0, i_mid, j_mid, configData.steps[1]);
            calcTempDiffRegion(model, configData, tempGrad, 0, i_mid, j_mid, configData.steps[1]);
        });
        std::thread th3([&](){
            calcGrainDiffEnergyRegion(model, configData, grainDiffEn, i_mid, configData.steps[0], 0, j_mid);
            calcPhaseDiffEnergyRegion(model, configData, phaseDiffEn, i_mid, configData.steps[0], 0, j_mid);
            calcTempDiffRegion(model, configData, tempGrad, i_mid, configData.steps[0], 0, j_mid);
        });
        std::thread th4([&](){
            calcGrainDiffEnergyRegion(model, configData, grainDiffEn, i_mid, configData.steps[0], j_mid, configData.steps[1]);
            calcPhaseDiffEnergyRegion(model, configData, phaseDiffEn, i_mid, configData.steps[0], j_mid, configData.steps[1]);
            calcTempDiffRegion(model, configData, tempGrad, i_mid, configData.steps[0], j_mid, configData.steps[1]);
        });

        th1.join();
        th2.join();
        th3.join();
        th4.join();

        // Now update the entire grid in a single call
        model.update(phaseDiffEn, tempGrad, grainDiffEn, configData, 0, configData.steps[0], 0, configData.steps[1]);
        model.resizeGrainDiffEn(grainDiffEn, model.numGrains);
    }

    std::string fileNameTemp = "TempGrid";
    std::ofstream output_file1(fileNameTemp);
    if (!output_file1.is_open()) {
        std::cerr << "Error: Could not open file " << fileNameTemp << std::endl;
    }
    for (int i=0;i < configData.steps[0];i++) {
        for (int j=0;j < configData.steps[1];j++) {
            output_file1 << model.grid[i][j].temp << ',';
        }
        output_file1 << std::endl;
    }
    output_file1.close();
    std::string fileNamePhase = "PhaseGrid";
    std::ofstream output_file2(fileNamePhase);
    if (!output_file2.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePhase << std::endl;
    }
    for (int i=0;i < configData.steps[0];i++) {
        for (int j=0;j < configData.steps[1];j++) {
            output_file2 << model.grid[i][j].phase << ',';
        }
        output_file2  << std::endl;
    }
    std::cout << "Calc Completed, Saved Data";
    output_file2.close();

    std::string fileNameGrain = "GrainGrid";
    std::ofstream output_file3(fileNameGrain);
    for (int g = 0; g < model.numGrains;g++) {
        for (int i=0;i < configData.steps[0];i++) {
            for (int j=0;j < configData.steps[1];j++) { 
                output_file3 << model.grid[i][j].grainPhases[g] << ',';
            }
            output_file3 << std::endl;
        }
        output_file3 << std::endl << std::endl;
    }
    output_file3.close();

    return 0;
}

// Helper function to resize grainDiffEn everywhere
