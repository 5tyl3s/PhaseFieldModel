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
    std::vector<int> nucLoc = {2,2};

    model.addGrain(nucLoc,configData);



    for (int t = 0; t < 200; t++) {
        grainDiffEn = calcGrainDiffEnergy(model,configData);
        phaseDiffEn = calcPhaseDiffEnergy(model,configData);
        tempGrad = calcTempDiff(model,configData);

        model.update(phaseDiffEn,tempGrad,grainDiffEn,configData); 
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