#pragma once
#include <vector>
struct fields {
    std::vector<int> steps; //[width,height]
    std::vector<std::vector<double>> grains;
    std::vector<std::vector<double>> orientation; //Euler Angles
    std::vector<double> temp;
    std::vector<double> phase;
    std::vector<double> particleComp;
};

struct config {
    double dx; //um
    double dt; //s
    double kSolid; //W/m*k
    double kLiquid = kSolid * 0.45; //W/m*k
    int success;
    double startTemp; //K
    double cellWidth;//um
    double cellHeight;//um
    fields modelFields;

};
config initModel();