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
    std::vector<std::vector<std::vector<double>>> grainDiffEn;
    std::vector<std::vector<double>> phaseDiffEn;
    std::vector<std::vector<double>> tempGrad;


    for (int t = 0; t < 5; t++) {
        std::cout << t;
        grainDiffEn = calcGrainDiffEnergy(model,configData);
        std::cout << "GrainEnDone\n";

        phaseDiffEn = calcPhaseDiffEnergy(model,configData);
        std::cout << "PhaseEnDone\n";
        tempGrad = calcTempDiff(model,configData);
        std::cout << "TempDone\n";
    
    }
    std::string fileName = "TempGrid";
    std::ofstream output_file(fileName);
    if (!output_file.is_open()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
    }
    for (int i=0;i < configData.steps[0]-1;i++) {
        for (int j=0;j < configData.steps[1]-2;j++) {
            output_file << model.grid[i][j].temp << ',';
        }
        output_file << model.grid[i][configData.steps[1]-1].temp << std::endl;
    }
    std::cout << "Calc Completed, Saved Data";
 
 


   
    return 0;

}