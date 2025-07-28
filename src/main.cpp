#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"

int main() {
    std::cout << "PhaseField Sim Starting" << std::endl;
    config configData = initModel();
    if (configData.success == 0) {
        std::cerr << "Initialization Failed";
        return 1;
    }
   
    return 0;

}