#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"

#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>

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
    auto start = std::chrono::steady_clock::now();
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

    // Calculate section sizes for 4x4 grid
    int i_eighth = configData.steps[0] / 8;
    int j_eighth = configData.steps[1] / 8;

    for (int t = 0; t < configData.timeSteps; t++) {
        if (t%100 == 0) std::cout << t << "/" << configData.timeSteps << std::endl;
        std::vector<std::thread> threads;

        // Create 16 threads in a 4x4 grid
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                int i_start = i * i_eighth;
                int i_end = (i == 7) ? configData.steps[0] : (i + 1) * i_eighth;
                int j_start = j * j_eighth;
                int j_end = (j == 7) ? configData.steps[1] : (j + 1) * j_eighth;

                // Launch phase and temp calculations
                threads.push_back(std::thread([&, i_start, i_end, j_start, j_end]() {
                    calcPhaseDiffEnergyRegion(model, configData, phaseDiffEn, i_start, i_end, j_start, j_end);
                }));


                // Launch grain calculations
                threads.push_back(std::thread([&, i_start, i_end, j_start, j_end]() {
                    calcGrainDiffEnergyRegion(model, configData, grainDiffEn, i_start, i_end, j_start, j_end);
                }));
            }
        }
        // Launch temperature calculation in main thread
        threads.push_back(std::thread([&]() {
            tempGrad = updateTemp(configData.tGrad, configData.coolingRate, t, configData.dx, configData.dt, configData.steps[0], configData.steps[1], configData.startTemp,  configData.minTemp);
        }));

        // Join all threads
        for (auto& thread : threads) {
            thread.join();
        }
        //std::cout << "Threads Joined" << std::endl;

        // Update the entire grid in a single call
        model.update(phaseDiffEn, tempGrad, grainDiffEn, configData, 0, configData.steps[0], 0, configData.steps[1]);
       // std::cout<<"Grid Updated"<<std::endl;
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
    double tempOut = 0;
    int gNum;
    
    for (int i=0;i < configData.steps[0];i++) {
        for (int j=0;j < configData.steps[1];j++) { 
            gNum = 0;
            tempOut = 0;
            for (int g = 0; g < model.grid[i][j].activeGrains.size(); g++) {
                if (model.grid[i][j].grainPhases[model.grid[i][j].activeGrains[g]] > tempOut) {
                    tempOut = model.grid[i][j].grainPhases[model.grid[i][j].activeGrains[g]];
                    gNum = model.grid[i][j].activeGrains[g];
                }
                
            }
            output_file3 << gNum << ',';
        }
        output_file3 << std::endl;        
    }

    output_file3.close();

    
    std::cout << "Elapsed(ms)=" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() << std::endl;
    return 0;
}

// Helper function to resize grainDiffEn everywhere
