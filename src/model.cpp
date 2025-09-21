#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>





double sumOtherGrainsSquared(int notIndex,node Node, int numGrains) {
    double notThisGrain = 0;
    for (int gr = 0; gr < 9;gr++) {
        notThisGrain = notThisGrain + (Node.grainPhases[gr]*Node.grainPhases[gr]);
    }
    return notThisGrain-(Node.grainPhases[notIndex]*Node.grainPhases[notIndex]);
    
}






double ifLiq(double temp, double meltTemp) {
    //0 Above melt, 1 below
    return 0.5*(1-(tanh(1000000*((temp/meltTemp)-1))));


}

double calcGrainBoundaryEnergy(eulerAngles orient, std::vector<double> gradient) {
    std::array<double,3> gradNormal = {-1*gradient[1],1*gradient[0],0};

    std::array<double,3> surfPlaneVec = eulerRotate(orient,gradNormal);
    
    double angle110 = dotAngle(surfPlaneVec,std::array<double,3> {1,1,0});
    double angle111 = dotAngle(surfPlaneVec,std::array<double,3> {1,1,1});
    double angle100 = dotAngle(surfPlaneVec,std::array<double,3> {1,0,0});
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
    std::array<double,3> eulerRotate(eulerAngles orient, std::array<double,3> rotatedVector) {
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

double dotAngle(std::array<double,3> vec1, std::array<double,3> vec2) {
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

double calcPhaseDiffEnergy(node* nd, config mConfig) {
    double pha = nd->phase;
    double part = nd->particleComp;
    double liq = ifLiq(nd->temp, mConfig.meltTemp);
    double eLoc = mConfig.phaseCoefficient * (
        -2 * (1 - pha) * liq + 2 * pha * (1 - liq) 
    );
    double grainSum = 0.0;
    for (int g = 0; g < 9; g++) {
        grainSum += nd->exists*nd->grainPhases[nd->activeGrains[g]] * nd->grainPhases[nd->activeGrains[g]];
    }
    eLoc += mConfig.grainPreCo * (-2 * (1 - pha) * grainSum) + (2*pha*part*part*mConfig.particleSolidIntEnergy+(-2+2*pha)*(part)*(part)*mConfig.particleLiquidIntEnergy);
    double eGrad = (nd->neighbors[0]->phase +
             nd->neighbors[1]->phase +
             nd->neighbors[2]->phase +
             nd->neighbors[3]->phase - (nd->neighbors[0]->exists + nd->neighbors[1]->exists + nd->neighbors[2]->exists + nd->neighbors[3]->exists) * pha) * mConfig.phaseGradCo;
    double diffFree = eLoc + eGrad + nd->particleComp*nd->particleComp*2*nd->phase*mConfig.particleSolidIntEnergy; ;
    std::cout << "Phase Diff Energy: " << diffFree << std::endl;
    return diffFree;

}



std::array<double,9> calcGrainDiffEnergy(node* nd, config mConfig) {
    std::array<double,9> diffFree;
    std::array<double,4> grainGradVals {0,0,0,0};
    for (int g = 0; g < 9; g++) {
        double gra = nd->grainPhases[nd->activeGrains[g]];

        for (int gg = 0; gg < 9; gg++) {
            if (nd->activeGrains[g] == nd->neighbors[0]->activeGrains[gg]) {
                grainGradVals[0] = nd->neighbors[0]->grainPhases[gg];
            }
            if (nd->activeGrains[g] == nd->neighbors[1]->activeGrains[gg]) {
                grainGradVals[1] = nd->neighbors[1]->grainPhases[gg];
            }
            if (nd->activeGrains[g] == nd->neighbors[2]->activeGrains[gg]) {
                grainGradVals[2] = nd->neighbors[2]->grainPhases[gg];
            }
            if (nd->activeGrains[g] == nd->neighbors[3]->activeGrains[gg]) {
                grainGradVals[3] = nd->neighbors[3]->grainPhases[gg];
            }
        }
        double grainGrad = (
            grainGradVals[0] + grainGradVals[1] + grainGradVals[2] + grainGradVals[3] -
            (nd->neighbors[0]->exists + nd->neighbors[1]->exists + nd->neighbors[2]->exists + nd->neighbors[3]->exists) * nd->grainPhases[g]) * -1 * mConfig.grainGradCo;

        double comp = sumOtherGrainsSquared(nd->activeGrains[g], *nd, 9);
        double solid = 1.0 - ifLiq(nd->temp, mConfig.meltTemp);

        double gbEnergy = calcGrainBoundaryEnergy(nd->orientations[g], {
            0.5*(grainGradVals[0] + grainGradVals[1])-nd->grainPhases[g],
            0.5*(grainGradVals[2] + grainGradVals[3])-nd->grainPhases[g]
        });
        grainGrad = grainGrad * gbEnergy;

        double grainEnergy = mConfig.grainPreCo * (gra * gra * gra * gra - gra + (1000 / mConfig.cellArea) * gra * comp * gbEnergy * mConfig.grainIntWidth + 2 * gra * solid * solid);

        diffFree[g] = mConfig.dt * (grainGrad + grainEnergy);
        
    }
    return diffFree;
}


double calcTemp(node* nd, config& modelConfig, int t) {
    return modelConfig.startTemp + (nd->heightPos*modelConfig.tGrad*modelConfig.dx) - (t*modelConfig.coolingRate*modelConfig.dt);
}

double calcParticleCompDiff(node* nd, config& modelConfig) {
    double part = nd->particleComp;
    double pha = nd->phase;
    double en = 2*part*pha*pha*modelConfig.particleSolidIntEnergy
    + 2*part*(1-pha)*(1-pha)*modelConfig.particleLiquidIntEnergy;
    return en;
}



