#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"

int main() {
    std::cout << "PhaseField Sim Starting" << std::endl;
    config configData = inputConfig();
    if (configData.success == 0) {
        std::cerr << "Config Failed";
        return 1;
    }

    gridField model;
    model.init(configData);
    eulerAngles orientTest = {3.14159/4,3.14159/4,0};
    std::array<std::array<double,3>,2> test = calcGrainBoundaryEnergy(orientTest,std::vector<double> {2,0,0});
    std::cout << test[0][0] << test[0][1]<< test[0][2] << std::endl << test[1][0] << test [1][1] << test[1][2] << std::endl;
   
    return 0;

}