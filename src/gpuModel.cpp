#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>


#include <vector>
std::vector<double> SH_coeffs; // Spherical harmonic coefficients
int Lmax = 8;                   // Msaximum degree (match MATLAB)
#include <fstream>
#include <string>


__global__ void tempCalc(double* tempGrid, double startTemp, double tGrad, double dx, double coolingRate, double dt, int gridCols, int timeStep, int vectorLength) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int index = idy * gridCols + idx;
    if (idx >= gridCols || idy >= vectorLength) return;

    double heightPos = idy; // assuming each row corresponds to a height position

    tempGrid[index] = startTemp + (heightPos * tGrad * dx) - tGrad * coolingRate * (dt * timeStep);
}


__global__ void phaseDiffEn(double* phaseDiffEnGrid, int gridCols, double dt, double dx, double phasePreCo, double* temps, double meltTemp, double* phases, double* particleCompGrid, double molarVolume, double* grainSum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int index = idy * gridCols + idx;
    if (idx >= gridCols || idy >= 1) return;

    double phase = phases[index];
    double part = particleCompGrid[index];
    double temp = temps[index];

    // Compute phase energy term here
    double phaseEnergy = (2*phase-2)*0.5*(1-(tanh(1000000*((temp/meltTemp)-1))))*((-0.0212 *temp + 61.3952)/molarVolume) + (2*phase*(1-0.5*(1-(tanh(1000000*((temp/meltTemp)-1)))))*((0.0349*temp - 101.0704)/molarVolume)) + phaseCoefficient*(10000*((2* phase-2) * grainSum[index]+(2*phase* (1 - grainSum[index])))) + ((phase) * pow(part,2)*particleSolidIntEnergy*(0.66*4*3.14149 * pow(particleRadius,2)/((4/3)*3.14159*pow(particleRadius,3)))) + ((2*phase-2) * (pow(part,2))*particleLiquidIntEnergy*(0.66*4*3.14149 * pow(particleRadius,2)/((4/3)*3.14159*pow(particleRadius,3))));
    double lapEnergy = 0.0; // Placeholder for gradient energy term
    
    phaseDiffEnGrid[index] = (dt / (dx * dx)) * phaseEnergy * 1.6e-16;
}
__global__ void grainDiffEn(double* grainDiffEnGrid, int gridCols, int vectorLength, double grainGradCo, double dx, double meltTemp, double grainPreCo, double grainIntWidth, double* temps, int maxGrainsPerNode, int* grainsHereGrid, double* orientationsGrid) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int index = idy * gridCols + idx;
    if (idx >= gridCols || idy >= vectorLength) return;

    // For simplicity, assume each node has a fixed number of grains (maxGrainsPerNode)
    for (int g = 0; g < maxGrainsPerNode; g++) {
        // Compute grain energy terms here
        // This is a placeholder; actual implementation would involve neighbor access and calculations
        double gra = 0.0; // Placeholder for grain phase fraction
        double temp = temps[index];

        double grainGrad = 0.0; // Placeholder for grain gradient term
        double comp = 0.0;      // Placeholder for composition interaction term

        double notUC = 1.0 - underCool(temp, meltTemp);
        double grainEnergy = grainPreCo * ((pow(gra, 3) - gra) * underCool(temp, meltTemp) + (gra * pow(comp, 2) * grainIntWidth));
        grainGrad = safeClamp(grainGrad, -1e12, 1e12);
        grainEnergy = safeClamp(grainEnergy, -1e12, 1e12);
        grainDiffEnGrid[index * maxGrainsPerNode + g] = grainGrad + grainEnergy;
    }
}