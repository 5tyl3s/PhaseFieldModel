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
   
    return 0;

}