#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>






std::vector<std::vector<double>> calcPhaseDiffEnergy(gridField field, config modelConfig) {
    std::vector<std::vector<double>> diffFree(modelConfig.steps[0],std::vector<double>(modelConfig.steps[1]));
    double pha;
    double eLoc;
    double eGrad;
    for (int ii = 0; ii < modelConfig.steps[0]; ii++) {
        for (int jj = 0; jj < modelConfig.steps[1]; jj++) {
            eLoc = 0.01*modelConfig.phaseCoefficient*((-1*(pha-2)*(ifLiq(field.grid[ii][jj].temp,modelConfig.meltTemp))));
            eGrad = 0.005*((pha*(1-ifLiq(field.grid[ii][jj].temp,modelConfig.meltTemp)))+(2*pha*(modelConfig.particleSlowingCoefficient*field.grid[ii][jj].particleComp))+(field.grid[ii][jj].neighbors[0]->phase +field.grid[ii][jj].neighbors[1]->phase +field.grid[ii][jj].neighbors[2]->phase +field.grid[ii][jj].neighbors[3]->phase -4))*modelConfig.phaseGradCo;
            //pha = field.grid[ii][jj].phase;
            //std::cout << pha << "   ";

            
            diffFree[ii][jj] =eLoc + eGrad;
            //std::cout << "Gradient Energy: " << eGrad << " Local Energy: " << eLoc << " Total Energy: " << diffFree[ii][jj] << std::endl;
            if (diffFree[ii][jj] < 0) {
                diffFree[ii][jj] = 0;
            }
            

        }
 
    }
    return diffFree; 



}

std::vector<std::vector<std::vector<double>>> calcGrainDiffEnergy(gridField field, config modelConfig) {
    std::vector<std::vector<std::vector<double>>> diffFree(modelConfig.steps[0],std::vector<std::vector<double>>(modelConfig.steps[1],std::vector<double>(field.numGrains)));
    double gra;
    double graEnergy;
    for (int i = 0; i < modelConfig.steps[0]; i++) {
        for (int j = 0; j < modelConfig.steps[1]; j++) {
            for (int g = 0; g < field.numGrains; g++){
                //std::cout << "Calculating Grain Energy for Grain: " << g << " At Location: " << i << "," << j << std::endl;
                gra = field.grid[i][j].grainPhases[g];
                //std::cout << gra << std::endl;
                graEnergy = modelConfig.grainPreCo*((field.grid[i][j].grainPhases[g]*field.grid[i][j].grainPhases[g]*field.grid[i][j].grainPhases[g]-field.grid[i][j].grainPhases[g])+(2*gra*compOtherGrains(g,field.grid[i][j],field.numGrains))+(2*field.grid[i][j].grainPhases[g]*field.grid[i][j].particleComp*modelConfig.particleSlowingCoefficient))+(modelConfig.grainIntWidth*calcGrainBoundaryEnergy(field.orientations[g],{field.grid[i][j].neighbors[0]->phase+field.grid[i][j].neighbors[1]->phase - 2*field.grid[i][j].phase,field.grid[i][j].neighbors[2]->phase+field.grid[i][j].neighbors[3]->phase - 2*field.grid[i][j].phase}))*(field.grid[i][j].neighbors[0]->grainPhases[g] +field.grid[i][j].neighbors[1]->grainPhases[g] +field.grid[i][j].neighbors[2]->grainPhases[g] +field.grid[i][j].neighbors[3]->grainPhases[g] -4);
                diffFree[i][j][g] = diffFree[i][j][g]+graEnergy;
                //std::cout << diffFree[i][j][g] << std::endl;
            }

        }
    }
    return diffFree;
}


std::vector<std::vector<double>> calcTempDiff(gridField& field, const config& modelConfig) {
    std::vector<std::vector<double>> tempDiff(modelConfig.steps[0],std::vector<double>(modelConfig.steps[1]));
    for (int j = 0; j<modelConfig.steps[1];j++) {
        tempDiff[0][j] = (1/(modelConfig.dx*modelConfig.dx))*(field.grid[0][j].neighbors[0]->temp + field.grid[0][j].neighbors[1]->temp + field.grid[0][j].neighbors[3]->temp - (3*field.grid[0][j].temp));
    }


    for (int i = 1; i < modelConfig.steps[0]-1; i++) {
        for (int j = 0; j < modelConfig.steps[1]; j++) {
            tempDiff[i][j] =(1/(modelConfig.dx*modelConfig.dx))*(field.grid[i][j].neighbors[0]->temp+field.grid[i][j].neighbors[1]->temp+field.grid[i][j].neighbors[2]->temp+field.grid[i][j].neighbors[3]->temp-4*field.grid[i][j].temp);
        }
    }
    int bottomMostStep = modelConfig.steps[0]-1;

    for (int jjjj = 0; jjjj<modelConfig.steps[1];jjjj++) {
        //std::cout << "Plus: " << field.grid[bottomMostStep][jjjj].neighbors[0]->temp + field.grid[bottomMostStep][jjjj].neighbors[1]->temp + field.grid[bottomMostStep][jjjj].neighbors[2]->temp + modelConfig.basePlateTemp << " Minus: " << 4*field.grid[bottomMostStep][jjjj].temp << std::endl;
        
        tempDiff[bottomMostStep][jjjj] = (1/(modelConfig.dx*modelConfig.dx))*(field.grid[bottomMostStep][jjjj].neighbors[0]->temp + field.grid[bottomMostStep][jjjj].neighbors[1]->temp + field.grid[bottomMostStep][jjjj].neighbors[2]->temp + modelConfig.basePlateTemp - (4*field.grid[bottomMostStep][jjjj].temp));
        
    }
    //std::cout << tempDiff[bottomMostStep][2] << std::endl;
    return tempDiff;

}























double compOtherGrains(int notIndex,node Node, int numGrains) {
    double notThisGrain = 0;
    for (int gr = 0; gr < numGrains;gr++) {
        notThisGrain = notThisGrain + Node.grainPhases[gr];
    }
    return notThisGrain-Node.grainPhases[notIndex];
    
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
    return freeEnVec;


}

double calcLocalFree(node nodeLoc, config mConfig) {
    return mConfig.phasePreCo*(((1-nodeLoc.phase)*(1-nodeLoc.phase)*(ifLiq(nodeLoc.temp,mConfig.meltTemp)))+(nodeLoc.phase*nodeLoc.phase*(1-ifLiq(nodeLoc.temp,mConfig.meltTemp))));
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
    std::cerr << "Grain Boundary Calculation Failed See model.cpp";
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


