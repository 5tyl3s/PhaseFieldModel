#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>

void say_hello() {
    std::cout << "Hello From Model"<< std::endl;
}




std::vector<std::vector<double>> calcPhaseDiffEnergy(gridField field, config modelConfig) {
    std::vector<std::vector<double>> diffFree;
    double pha;
    for (int i = 0; i < modelConfig.steps[0]-1; i++) {
        for (int j = 0; j < modelConfig.steps[1]-1; j++) {
            pha = field.grid[i][j].phase;
            diffFree[i][j] = modelConfig.phaseCoefficient*(((2*pha-2)*(ifLiq(field.grid[i][j].temp,modelConfig.meltTemp)))+(2*pha*(1-ifLiq(field.grid[i][j].temp,modelConfig.meltTemp)))+(2*pha*(modelConfig.particleSlowingCoefficient*field.grid[i][j].particleComp)))+(field.grid[i][j].neighbors[0]->phase +field.grid[i][j].neighbors[1]->phase +field.grid[i][j].neighbors[2]->phase +field.grid[i][j].neighbors[3]->phase -4)*modelConfig.phaseGradCo;
            

        }
    }



}

std::vector<std::vector<std::vector<double>>> calcGrainDiffEnergy(gridField field, config modelConfig) {
    std::vector<std::vector<double>> diffFree;
    double gra;
    double graEnergy;
    for (int i = 0; i < modelConfig.steps[0]-1; i++) {
        for (int j = 0; j < modelConfig.steps[1]-1; j++) {
            for (int g = 0; g < field.numGrains; g++){
                gra = field.grid[i][j].grainPhases[g];
                graEnergy = modelConfig.grainPreCo*((field.grid[i][j].grainPhases[g]*field.grid[i][j].grainPhases[g]*field.grid[i][j].grainPhases[g]-field.grid[i][j].grainPhases[g])+(2*gra*compOtherGrains(g,field.grid[i][j]))+(2*field.grid[i][j].grainPhases[g]*field.grid[i][j].particleComp*modelConfig.particleSlowingCoefficient))+(modelConfig.grainIntWidth*calcGrainBoundaryEnergy(field.orientations[g],{field.grid[i][j].neighbors[0]->phase+field.grid[i][j].neighbors[1]->phase - 2*field.grid[i][j].phase,field.grid[i][j].neighbors[2]->phase+field.grid[i][j].neighbors[3]->phase - 2*field.grid[i][j].phase}))*(field.grid[i][j].neighbors[0]->grainPhases[g] +field.grid[i][j].neighbors[1]->grainPhases[g] +field.grid[i][j].neighbors[2]->grainPhases[g] +field.grid[i][j].neighbors[3]->grainPhases[g] -4);

                diffFree[i][j] = diffFree[i][j]+graEnergy;
            }

        }
    }
}

double compOtherGrains(int notIndex,node Node) {
    double notThisGrain = 0;
    for (int gr = 0; gr < size(Node.grainPhases);gr++) {
        notThisGrain = notThisGrain + Node.grainPhases[gr];
    }
    notThisGrain = notThisGrain-Node.grainPhases[notIndex];
}





//Free Energy Calcs, used for validating energy decrease??
std::vector<std::vector<double>> calcFreeEnergy(gridField field,config modelConfig) {
    std::vector<std::vector<double>> freeEnVec;

    for (int i = 0; i < modelConfig.steps[0]-1; i++) {
        for (int j = 0; j < modelConfig.steps[1]-1; j++) {
            double fLocal = calcLocalFree(field.grid[i][j],modelConfig);
            double fPhaseGrad = 0.5*modelConfig.phaseCoefficient*(field.grid[i][j].neighbors[0]->phase +field.grid[i][j].neighbors[1]->phase +field.grid[i][j].neighbors[2]->phase +field.grid[i][j].neighbors[3]->phase -4);//Von Nueman laplacian
            double gbEn = 0;
            for (int gr = 0; gr < field.numGrains; gr++)
                if (field.grid[i][j].grainPhases[gr] !=0) {
                    gbEn = gbEn + calcLocalFree(field.grid[i][j], modelConfig);
                };

            double fGrainGrad = gbEn*modelConfig.grainIntWidth;

            freeEnVec[i][j] = fLocal+fPhaseGrad+fGrainGrad;


        }
    }


}

double calcLocalFree(node nodeLoc, config mConfig) {
    double energyPhase = mConfig.phasePreCo*(((1-nodeLoc.phase)*(1-nodeLoc.phase)*(ifLiq(nodeLoc.temp,mConfig.meltTemp)))+(nodeLoc.phase*nodeLoc.phase*(1-ifLiq(nodeLoc.temp,mConfig.meltTemp))));


}
double ifLiq(double temp, double meltTemp) {
    //0 Above melt, 1 below
    double liqVal = 0.5*(1-(tanh((temp/meltTemp)-1)));
    return liqVal;

}

double calcGrainBoundaryEnergy(eulerAngles orient, std::array<double,2> gradient) {
    std::array<double, 3> gradNormal = {-1*gradient[1],1*gradient[0],0};
    std::array<double, 3> surfPlaneVec =eulerRotate(orient,gradNormal);
    std::cout << surfPlaneVec[0] << ' ' << surfPlaneVec[1] << ' ' << surfPlaneVec[2] << std::endl;

    double angle110 = dotAngle(surfPlaneVec,std::array<double,3> {1,1,0});
    double angle111 = dotAngle(surfPlaneVec,std::array<double,3> {1,1,1});
    double angle100 = dotAngle(surfPlaneVec,std::array<double,3> {1,0,0});

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
    double adotb = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
    double mag1 = std::sqrt((vec1[0]*vec1[0])+(vec1[1]*vec1[1])+(vec1[2]*vec1[2]));
    double mag2 = std::sqrt((vec2[0]*vec2[0])+(vec2[1]*vec2[1])+(vec2[2]*vec2[2]));
    double angle = acos(adotb/(mag1*mag2));
    while (angle > 90) {
        angle = angle - 90;
    }
    while (angle < 0) {
        angle = angle + 90;
    }
    return angle;

}