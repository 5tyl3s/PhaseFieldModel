#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

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
    std::vector<std::vector<double>> tempGrad;
    std::vector<int> nucLoc = {37,20};

    


    
    for (int t = 0; t < 10; t++) {
        std::cout << t << std::endl;
        grainDiffEn = calcGrainDiffEnergy(model,configData);
        phaseDiffEn = calcPhaseDiffEnergy(model,configData);
        tempGrad = calcTempDiff(model,configData);

        model.update(phaseDiffEn,tempGrad,grainDiffEn,configData); 
    }
    model.addGrain(nucLoc,configData);
    for (int t = 0; t < 500; t++) {
        std::cout << t << std::endl;
        grainDiffEn = calcGrainDiffEnergy(model,configData);
        phaseDiffEn = calcPhaseDiffEnergy(model,configData);
        tempGrad = calcTempDiff(model,configData);

        model.update(phaseDiffEn,tempGrad,grainDiffEn,configData); 
    }








    std::string fileNameTemp = "TempGrid";
    std::ofstream output_file1(fileNameTemp);
    if (!output_file1.is_open()) {
        std::cerr << "Error: Could not open file " << fileNameTemp << std::endl;
    }
    for (int i=0;i < configData.steps[0]-1;i++) {
        for (int j=0;j < configData.steps[1]-2;j++) {
            output_file1 << model.grid[i][j].temp << ',';
        }
        output_file1 << model.grid[i][configData.steps[1]-1].temp << std::endl;
    }
    output_file1.close();
    std::string fileNamePhase = "PhaseGrid";
    std::ofstream output_file2(fileNamePhase);
    if (!output_file2.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePhase << std::endl;
    }
    for (int i=0;i < configData.steps[0]-1;i++) {
        for (int j=0;j < configData.steps[1]-2;j++) {
            output_file2 << model.grid[i][j].phase << ',';
        }
        output_file2 << model.grid[i][configData.steps[1]-1].phase << std::endl;
    }
    std::cout << "Calc Completed, Saved Data";
    output_file2.close();

    std::string fileNameGrain = "GrainGrid";
    std::ofstream output_file3(fileNameGrain);
    for (int g = 0; g < model.numGrains;g++) {
        for (int i=0;i < configData.steps[0]-1;i++) {
            for (int j=0;j < configData.steps[1]-2;j++) { 
                
                output_file3 << model.grid[i][j].grainPhases[g] << ',';

            }
            output_file3 << std::endl;
                
        }
        output_file3 << std::endl << std::endl;
    }
    output_file3.close();
 
 


   
    return 0;

}