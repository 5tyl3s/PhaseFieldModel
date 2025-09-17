#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include "importConfig.hpp"
#include <omp.h>
#include <direct.h> // For _getcwd on Windows

#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>

// These functions must be implemented to fill only the specified region of the output array
int main() {
    std::cout << "Program started" << std::endl << std::flush;
    
    // Print working directory
    char cwd[256];
    if (_getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cout << "Current working directory: " << cwd << std::endl << std::flush;
    } else {
        std::cerr << "Could not get current working directory" << std::endl << std::flush;
    }

    std::cout << "About to start PhaseField Sim..." << std::endl << std::flush;
    auto start = std::chrono::steady_clock::now();
    
    std::cout << "Attempting to read config file..." << std::endl << std::flush;
    config configData = inputConfig("config/config.txt");
    if (configData.success == 0) {
        std::cerr << "Config Failed to load from: config/config.txt" << std::endl << std::flush;
        return 1;
    }
    std::cout << "Config loaded successfully" << std::endl << std::flush;



    gridField model;
    model.init(configData);
    std::array<std::array<double,9>,1000000> grainDiffEn;
    std::array<double,1000000> phaseDiffEn;
    std::array<double,1000000>  tempGrad;
    std::array<double,1000000> tempPartComp;
    node* nd;



    for (int t = 0; t < configData.timeSteps; t++) {
        if (t%100 == 0) std::cout << t << "/" << configData.timeSteps << std::endl;
        #pragma omp parallel for
        for (int node = 0; node < 1000000; node++) {
            nd = model.allNodes[node];

    
            phaseDiffEn[node] = calcPhaseDiffEnergy(nd, configData);
            tempPartComp[node] = calcParticleCompDiff(nd, configData);
            grainDiffEn[node] = calcGrainDiffEnergy(nd, configData);
        }
        tempGrad = updateTemp(configData.tGrad, configData.coolingRate, t, configData.dx, configData.dt, configData.steps[0], configData.steps[1], configData.startTemp, configData.minTemp);

        


        // Join all threads

        //std::cout << "Threads Joined" << std::endl;

        // Update the entire grid in a single call
        model.update(phaseDiffEn, tempGrad, grainDiffEn, tempPartComp);

    }

    std::string fileNameTemp = "TempGrid.csv";
    std::ofstream output_file1(fileNameTemp);
    if (!output_file1.is_open()) {
        std::cerr << "Error: Could not open file " << fileNameTemp << std::endl;
    }
    for (int i=0;i < configData.steps[0];i++) {
        for (int j=0;j < configData.steps[1];j++) {
            output_file1 << model.grid[i][j].temp;
            if (j < configData.steps[1] - 1) {
                output_file1 << ',';
            }
        }
        output_file1 << std::endl;
    }
    output_file1.close();
    std::string fileNamePhase = "PhaseGrid.csv";
    std::ofstream output_file2(fileNamePhase);
    if (!output_file2.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePhase << std::endl;
    }
    for (int i=0;i < configData.steps[0];i++) {
        for (int j=0;j < configData.steps[1];j++) {
            output_file2  << model.grid[i][j].phase;
            if (j < configData.steps[1] - 1) {
                output_file2  << ',';
            }
        }
        output_file2  << std::endl;
    }
    std::cout << "Calc Completed, Saved Data";
    output_file2.close();

    std::string fileNameGrain = "GrainGrid.csv";
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
                    gNum = model.grid[i][j].activeGrains[g]+1;
                }
                
            }
            output_file3 << gNum;
            if (j < configData.steps[1] - 1) {
                output_file3 << ',';
            }
            
        }
        output_file3 << std::endl;        
    }

    output_file3.close();


    std::string fileNamePart = "ParticleGrid.csv";
    std::ofstream output_file4(fileNamePart);
    if (!output_file4.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePart << std::endl;
    }
    for (int i=0;i < configData.steps[0];i++) {
        for (int j=0;j < configData.steps[1];j++) {
            output_file4  << model.grid[i][j].particleComp;
            if (j < configData.steps[1] - 1) {
                output_file4  << ',';
            }
        }
        output_file4  << std::endl;
    }
    std::cout << "Calc Completed, Saved Data";
    output_file4   .close();

    
    std::cout << "Elapsed(ms)=" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() << std::endl;
    return 0;

}

// Helper function to resize grainDiffEn everywhere

