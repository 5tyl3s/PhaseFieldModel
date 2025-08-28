#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>





double sumOtherGrainsSquared(int notIndex,node Node, int numGrains) {
    double notThisGrain = 0;
    for (int gr = 0; gr < numGrains;gr++) {
        notThisGrain = notThisGrain + (Node.grainPhases[gr]*Node.grainPhases[gr]);
    }
    return notThisGrain-(Node.grainPhases[notIndex]*Node.grainPhases[notIndex]);
    
}






double ifLiq(double temp, double meltTemp) {
    //0 Above melt, 1 below
    return 0.5*(1-(tanh(1000000*((temp/meltTemp)-1))));


}

double calcGrainBoundaryEnergy(eulerAngles orient, std::vector<double> gradient) {
    std::vector<double> gradNormal = {-1*gradient[1],1*gradient[0],0};

    std::vector<double> surfPlaneVec =eulerRotate(orient,gradNormal);
    
    double angle110 = dotAngle(surfPlaneVec,std::vector<double> {1,1,0});
    double angle111 = dotAngle(surfPlaneVec,std::vector<double> {1,1,1});
    double angle100 = dotAngle(surfPlaneVec,std::vector<double> {1,0,0});
    if (std::isnan(angle110)) { 
        return 0;
    }


    if (angle110 < angle111) {
        if (angle110 < angle100) {
            return 2.783;
        }
    }
    if (angle111<angle100) {
        if (angle111<angle110) {
            return 2.962;
        }
    } else {
        return 3.192;
    }
    std::cerr << gradNormal[0] << " " << gradNormal[1] << " " << gradNormal[2] << std::endl;

    std::cerr << surfPlaneVec[0] << surfPlaneVec[1] << surfPlaneVec[2] << std::endl;
    std::cerr << angle110 << std::endl;
    std::cerr << angle111 << std::endl;
    std::cerr << angle100 << std::endl;
    std::cerr << "Grain Boundary Calculation Failed See model.cpp" << std::endl;
    return 999;



}














//Math Equations
std::vector<double> eulerRotate(eulerAngles orient, std::vector<double> rotatedVector) {
    double c1 = std::cos(orient.theta1), s1 = std::sin(orient.theta1);
    double cP = std::cos(orient.phi),  sP = std::sin(orient.phi);
    double c2 = std::cos(orient.theta2), s2 = std::sin(orient.theta2);

    // Build rotation matrix R = Rz(phi1) * Rx(Phi) * Rz(phi2)
    // (row-major order)
    double R00 =  c1*c2 - s1*s2*cP;
    double R01 = -c1*s2 - s1*c2*cP;
    double R02 =  s1*sP;

    double R10 =  s1*c2 + c1*s2*cP;
    double R11 = -s1*s2 + c1*c2*cP;
    double R12 = -c1*sP;

    double R20 =  s2*sP;
    double R21 =  c2*sP;
    double R22 =  cP;

    // Multiply R * v_in
    return {
        R00 * rotatedVector[0] + R01 * rotatedVector[1] + R02 * rotatedVector[2],
        R10 * rotatedVector[0] + R11 * rotatedVector[1] + R12 * rotatedVector[2],
        R20 * rotatedVector[0] + R21 * rotatedVector[1] + R22 * rotatedVector[2]
    };


}

double dotAngle(std::vector<double> vec1, std::vector<double> vec2) {
    //std::cout << "Vector 1: " << vec1[0] << ", " << vec1[1] << ", " << vec1[2] << std::endl;
    //std::cout << "Vector 2: " << vec2[0] << ", " << vec2[1] << ", " << vec2[2] << std::endl;
    double adotb = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
    //std::cout << "Dot Product: " << adotb << std::endl;
    double mag1 = std::sqrt((vec1[0]*vec1[0])+(vec1[1]*vec1[1])+(vec1[2]*vec1[2]));
    //std::cout << "Magnitude 1: " << mag1 << std::endl;
    double mag2 = std::sqrt((vec2[0]*vec2[0])+(vec2[1]*vec2[1])+(vec2[2]*vec2[2]));
    //std::cout << "Magnitude 2: " << mag2 << std::endl;
    //std::cout << "Acos of this: " << (adotb/(mag1*mag2)) << std::endl;
    double angle = acos(adotb/(mag1*mag2));

    while (angle > 90) {
        angle = angle - 90;
    }
    while (angle < 0) {
        angle = angle + 90;
    }
    //std::cout << angle << std::endl;
    return angle;

}

void calcPhaseDiffEnergyRegion(gridField& field, config& modelConfig,
    std::vector<std::vector<double>>& diffFree,
    int i_start, int i_end, int j_start, int j_end) {
    double pha, eLoc, eGrad;
    for (int ii = i_start; ii < i_end; ii++) {
        for (int jj = j_start; jj < j_end; jj++) {
            pha = field.grid[ii][jj].phase;
            double liq = ifLiq(field.grid[ii][jj].temp, modelConfig.meltTemp);
            eLoc = modelConfig.phaseCoefficient * (
                -2 * (1 - pha) * liq + 2 * pha * (1 - liq)
            );
            double grainSum = 0.0;
            for (int g = 0; g < field.grid[ii][jj].activeGrains.size(); g++) {
                grainSum += field.grid[ii][jj].grainPhases[field.grid[ii][jj].activeGrains[g]] * field.grid[ii][jj].grainPhases[field.grid[ii][jj].activeGrains[g]];
                //std::cout << "Grain " << field.grid[ii][jj].activeGrains[g] << " Phase: " << field.grid[ii][jj].grainPhases[field.grid[ii][jj].activeGrains[g]] << std::endl;
                //std::cout << "Grain Sum: " << grainSum << std::endl;
            }

            
            
            eLoc += modelConfig.grainPreCo * (-2 * (1 - pha) * grainSum);
            eLoc += 2 * pha * (modelConfig.particleSlowingCoefficient * field.grid[ii][jj].particleComp);
            eGrad = (field.grid[ii][jj].neighbors[0]->phase +
                     field.grid[ii][jj].neighbors[1]->phase +
                     field.grid[ii][jj].neighbors[2]->phase +
                     field.grid[ii][jj].neighbors[3]->phase - 4 * pha) * modelConfig.phaseGradCo;
            diffFree[ii][jj] = eLoc + eGrad;
            //std::cout << "Phase Energy at (" << ii << "," << jj << "): " << diffFree[ii][jj] << " From Grad: " << eGrad << std::endl;
        }
    }
}

void calcGrainDiffEnergyRegion(gridField& field, config& modelConfig,
    std::vector<std::vector<std::vector<double>>>& diffFree,
    int i_start, int i_end, int j_start, int j_end) {
    for (int i = i_start; i < i_end; i++) {
        for (int j = j_start; j < j_end; j++) {
            for (int g = 0; g < field.grid[i][j].activeGrains.size(); g++) {
                double gra = field.grid[i][j].grainPhases[field.grid[i][j].activeGrains[g]];
                double grainGrad = (
                    field.grid[i][j].neighbors[0]->grainPhases[field.grid[i][j].activeGrains[g]] +
                    field.grid[i][j].neighbors[1]->grainPhases[field.grid[i][j].activeGrains[g]] +
                    field.grid[i][j].neighbors[2]->grainPhases[field.grid[i][j].activeGrains[g]] +
                    field.grid[i][j].neighbors[3]->grainPhases[field.grid[i][j].activeGrains[g]] -
                    4 * gra
                ) * -1 * modelConfig.grainGradCo;

                double comp = sumOtherGrainsSquared(field.grid[i][j].activeGrains[g], field.grid[i][j], field.numGrains);
                double solid = 1.0 - ifLiq(field.grid[i][j].temp, modelConfig.meltTemp);

                double gbEnergy = calcGrainBoundaryEnergy(field.orientations[field.grid[i][j].activeGrains[g]], {
                    0.5 * field.grid[i][j].neighbors[0]->grainPhases[field.grid[i][j].activeGrains[g]] + 0.5 * field.grid[i][j].neighbors[1]->grainPhases[field.grid[i][j].activeGrains[g]] - field.grid[i][j].grainPhases[field.grid[i][j].activeGrains[g]],
                    0.5 * field.grid[i][j].neighbors[2]->grainPhases[field.grid[i][j].activeGrains[g]] + 0.5 * field.grid[i][j].neighbors[3]->grainPhases[field.grid[i][j].activeGrains[g]] - field.grid[i][j].grainPhases[field.grid[i][j].activeGrains[g]],
                });
                grainGrad = grainGrad * gbEnergy;

                double grainEnergy = modelConfig.grainPreCo * (gra * gra * gra * gra - gra + (1000 / modelConfig.cellArea) * gra * comp * gbEnergy * modelConfig.grainIntWidth + 2 * gra * solid * solid);

                diffFree[i][j][field.grid[i][j].activeGrains[g]] = modelConfig.dt * (grainGrad + grainEnergy);
            }
        }
    }
}

void calcTempDiffRegion(gridField& field, config& modelConfig,
    std::vector<std::vector<double>>& tempDiff,
    int i_start, int i_end, int j_start, int j_end) {
    // Top boundary
    if (i_start == 0) {
        for (int j = j_start; j < j_end; j++) {
            tempDiff[0][j] = (1 / (modelConfig.dx * modelConfig.dx)) *
                (field.grid[0][j].neighbors[0]->temp +
                 field.grid[0][j].neighbors[1]->temp +
                 field.grid[0][j].neighbors[3]->temp -
                 (3 * field.grid[0][j].temp));
        }
    }

    // Interior
    for (int i = std::max(i_start, 1); i < std::min(i_end, modelConfig.steps[0] - 1); i++) {
        for (int j = j_start; j < j_end; j++) {
            tempDiff[i][j] = (1 / (modelConfig.dx * modelConfig.dx)) *
                (field.grid[i][j].neighbors[0]->temp +
                 field.grid[i][j].neighbors[1]->temp +
                 field.grid[i][j].neighbors[2]->temp +
                 field.grid[i][j].neighbors[3]->temp -
                 4 * field.grid[i][j].temp);
        }
    }

    // Bottom boundary
    int bottomMostStep = modelConfig.steps[0] - 1;
    if (i_end == modelConfig.steps[0]) {
        for (int j = j_start; j < j_end; j++) {
            tempDiff[bottomMostStep][j] = (1 / (modelConfig.dx * modelConfig.dx)) *
                (field.grid[bottomMostStep][j].neighbors[0]->temp +
                 field.grid[bottomMostStep][j].neighbors[1]->temp +
                 field.grid[bottomMostStep][j].neighbors[2]->temp +
                 modelConfig.basePlateTemp -
                 (4 * field.grid[bottomMostStep][j].temp));
        }
    }
}




