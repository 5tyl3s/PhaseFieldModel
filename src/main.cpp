#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include "importConfig.hpp"
#include "gridField.hpp"
#include <omp.h>
#include <direct.h> // For _getcwd on Windows

#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>

// Use compile-time TOTAL_NODES from header so sizes match the grid definition
std::array<double, TOTAL_NODES> newTemp;
std::array<std::array<double,9>, TOTAL_NODES> grainDiffEn;
std::array<double, TOTAL_NODES> phaseDiffEn;
std::array<double, TOTAL_NODES> tempGrad;
std::array<double, TOTAL_NODES> tempPartComp;

int main() {
    std::cout << "Hello, World!" << std::endl;
    std::cout << "Attempting to read config file..." << std::endl << std::flush;
    config configData = inputConfig("config/config.txt");
    if (configData.success == 0) {
        std::cerr << "Config Failed to load from: config/config.txt" << std::endl << std::flush;
        return 1;
    }

    //gridField model2;
    globalField.init(configData);
    char cwd2[256];
    if (_getcwd(cwd2, sizeof(cwd2)) != NULL) {
        std::cout << "Current working directory: " << cwd2 << std::endl << std::flush;
    } else {
        std::cerr << "Could not get current working directory" << std::endl << std::flush;
    }    std::cout << "About to start PhaseField Sim..." << std::endl << std::flush;
    auto start1 = std::chrono::steady_clock::now();
    node* nd2;
    std::cout << "Starting Time Steps..." << std::endl << std::flush;

    for (int t = 0; t < configData.timeSteps; t++) {
        if (t%10 == 0) std::cout << t << "/" << configData.timeSteps << std::endl;
        #pragma omp parallel for private(nd2)
        for (int node = 0; node < TOTAL_NODES; node++) {
            nd2 = globalField.allNodes[node];
            tempGrad[node] = calcTemp(nd2,configData,t);
            phaseDiffEn[node] = calcPhaseDiffEnergy(nd2, configData);
            tempPartComp[node] = calcParticleCompDiff(nd2, configData);
            grainDiffEn[node] = calcGrainDiffEnergy(nd2, configData);
        }
        

        globalField.update(phaseDiffEn, grainDiffEn, tempPartComp, tempGrad);

    }


    // Write output files using GRID_ROWS / GRID_COLS
    std::string fileNameTemp = "TempGrid.csv";
    std::ofstream output_file1(fileNameTemp);
    if (!output_file1.is_open()) {
        std::cerr << "Error: Could not open file " << fileNameTemp << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            output_file1 << /* value */ globalField.grid[i][j].temp;
            if (j < GRID_COLS - 1) output_file1 << ",";
        }
        output_file1 << "\n";
    }
    output_file1.close();
    std::string fileNamePhase = "PhaseGrid.csv";
    std::ofstream output_file2(fileNamePhase);
    if (!output_file2.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePhase << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            output_file2  << globalField.grid[i][j].phase;
            if (j < GRID_COLS - 1) output_file2  << ",";
        }
        output_file2  << "\n";
    }
    std::cout << "Calc Completed, Saved Data";
    output_file2.close();

    std::string fileNameGrain = "GrainGrid.csv";
    std::ofstream output_file3(fileNameGrain);
    double tempOut = 0;
    int gNum;
    
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) { 
            gNum = 0;
            tempOut = 0;
            for (int g = 0; g < 9; g++) {
                if (globalField.grid[i][j].grainPhases[g] > tempOut) {
                    tempOut = globalField.grid[i][j].grainPhases[g];
                    gNum = globalField.grid[i][j].activeGrains[g]+1;
                }
                
            }
            output_file3 << gNum;
            if (j < GRID_COLS - 1) {
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
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            output_file4  << globalField.grid[i][j].particleComp;
            if (j < GRID_COLS - 1) {
                output_file4  << ',';
            }
        }
        output_file4  << std::endl;
    }
    std::cout << "Calc Completed, Saved Data";
    output_file4   .close();

    
    std::cout << "Elapsed(ms)=" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start1).count() << std::endl;
    return 0;

}

// Helper function to resize grainDiffEn everywhere

