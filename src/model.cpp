#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <atomic>

#include <vector>
std::vector<double> SH_coeffs; // Spherical harmonic coefficients
int Lmax = 8;                   // Maximum degree (match MATLAB)

constexpr int SH_LUT_THETA_SIZE = 3600;
constexpr int SH_LUT_PHI_SIZE   = 3600;
std::vector<float> SH_energyLUT;

#include <fstream>
#include <string>

struct EnergyProfiler {
    std::atomic<long long> phaseThermoNs{0};
    std::atomic<long long> phaseLaplacianNs{0};
    std::atomic<long long> phaseTotalNs{0};
    std::atomic<long long> phaseCalls{0};

    std::atomic<long long> grainNeighborNs{0};
    std::atomic<long long> grainCountNs{0};
    std::atomic<long long> grainGradientNs{0};
    std::atomic<long long> grainGbEnergyNs{0};
    std::atomic<long long> grainEnergyNs{0};
    std::atomic<long long> grainTotalNs{0};
    std::atomic<long long> grainCalls{0};

    std::atomic<long long> particleMuNs{0};
    std::atomic<long long> particleTotalNs{0};
    std::atomic<long long> particleCalls{0};
} energyProfiler;

auto safeClamp = [](double v, double low, double high){
    if (std::isnan(v) || std::isinf(v)) return (low+high)/2.0;
    return std::max(low, std::min(high, v));
};

double factorial(int n);
double P_lm(int l, int m, double x);
double Y_lm_real(int l, int m, double theta, double phi);

bool readCoeffs(const std::string&filename) {
    std::ifstream file(filename);
    std::cout << "Reading Spherical Harmonic Coefficients from: " << filename << std::endl;
    if(!file.is_open()) return false;

    SH_coeffs.clear();
    std::cout << "Loading Coefficients..." << std::endl;
    double val;
    while(file >> val) SH_coeffs.push_back(val);
    std::cout << "Loaded " << SH_coeffs.size() << " Coefficients." << std::endl;
    file.close();

    int expected_count = (Lmax + 1) * (Lmax + 1);
    if ((int)SH_coeffs.size() < expected_count) {
        std::cerr << "SH_coeffs size " << SH_coeffs.size() << " < expected " << expected_count << std::endl;
        return false;
    }

    SH_energyLUT.clear();
    SH_energyLUT.resize(SH_LUT_THETA_SIZE * SH_LUT_PHI_SIZE);
    const double PI = 3.141592653589793;
    const double dTheta = PI / SH_LUT_THETA_SIZE;
    const double dPhi = 2.0 * PI / SH_LUT_PHI_SIZE;

    for (int ti = 0; ti < SH_LUT_THETA_SIZE; ++ti) {
        double theta = dTheta * (ti + 0.5);
        for (int pi = 0; pi < SH_LUT_PHI_SIZE; ++pi) {
            double phi = dPhi * pi;
            double energy = 0.0;
            int idx = 0;
            for (int l = 0; l <= Lmax; ++l) {
                for (int m = -l; m <= l; ++m) {
                    energy += SH_coeffs[idx++] * Y_lm_real(l, m, theta, phi);
                }
            }
            if (!std::isfinite(energy)) energy = 2.9297;
            energy = std::max(2.9297, std::min(8.0, energy));
            SH_energyLUT[pi + ti * SH_LUT_PHI_SIZE] = static_cast<float>(energy);
        }
    }

    std::cout << "Built spherical harmonic energy lookup table (" << SH_LUT_THETA_SIZE << "x" << SH_LUT_PHI_SIZE << ")." << std::endl;
    return true;
}

// factorial helper
double factorial(int n) {
    double f = 1.0;
    for(int i=2;i<=n;i++) f *= i;
    return f;
}

// Associated Legendre Polynomial P_l^m(x)
double P_lm(int l, int m, double x) {
    if(m < 0) m = -m;
    double pmm = 1.0;
    if(m > 0) {
        double somx2 = std::sqrt(1-x*x);
        double fact = 1.0;
        for(int i=1;i<=m;i++) { pmm *= -fact*somx2; fact += 2; }
    }
    if(l == m) return pmm;
    double pmmp1 = x*(2*m+1)*pmm;
    if(l == m+1) return pmmp1;
    double pll = 0.0;
    for(int ll = m+2; ll <= l; ll++) {
        pll = ((2*ll-1)*x*pmmp1 - (ll+m-1)*pmm)/(ll-m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

// Real spherical harmonic Y_lm
double Y_lm_real(int l, int m, double theta, double phi) {
    double norm = std::sqrt((2*l+1)/(4*3.141592653589793) * factorial(l-m)/factorial(l+m));
    if(m==0) return norm * P_lm(l,0,std::cos(theta));
    if(m>0) return std::sqrt(2.0)*norm*P_lm(l,m,std::cos(theta))*std::cos(m*phi);
    else return std::sqrt(2.0)*norm*P_lm(l,-m,std::cos(theta))*std::sin(-m*phi);
}

double underCool(double temp, double meltTemp) {
    return 0.5*(1-(tanh(1000000*((temp/meltTemp)-1))));
}

void resetEnergyProfilingStats() {
    energyProfiler.phaseThermoNs = 0;
    energyProfiler.phaseLaplacianNs = 0;
    energyProfiler.phaseTotalNs = 0;
    energyProfiler.phaseCalls = 0;
    energyProfiler.grainNeighborNs = 0;
    energyProfiler.grainCountNs = 0;
    energyProfiler.grainGradientNs = 0;
    energyProfiler.grainGbEnergyNs = 0;
    energyProfiler.grainEnergyNs = 0;
    energyProfiler.grainTotalNs = 0;
    energyProfiler.grainCalls = 0;
    energyProfiler.particleMuNs = 0;
    energyProfiler.particleTotalNs = 0;
    energyProfiler.particleCalls = 0;
}

void printEnergyProfilingStats() {
    auto printMs = [](long long ns) {
        return static_cast<double>(ns) * 1e-6;
    };

    std::cout << "\n=== Energy Profiling Results ===" << std::endl;
    if (energyProfiler.phaseCalls > 0) {
        std::cout << "Phase differential calls: " << energyProfiler.phaseCalls << std::endl;
        std::cout << "  Thermodynamic energy time: " << printMs(energyProfiler.phaseThermoNs.load()) << " ms" << std::endl;
        std::cout << "  Laplacian time: " << printMs(energyProfiler.phaseLaplacianNs.load()) << " ms" << std::endl;
        std::cout << "  Total phase diff time: " << printMs(energyProfiler.phaseTotalNs.load()) << " ms" << std::endl;
    }
    if (energyProfiler.grainCalls > 0) {
        std::cout << "Grain differential calls: " << energyProfiler.grainCalls << std::endl;
        std::cout << "  Neighbor gather time: " << printMs(energyProfiler.grainNeighborNs.load()) << " ms" << std::endl;
        std::cout << "  Grain count/laplacian time: " << printMs(energyProfiler.grainCountNs.load()) << " ms" << std::endl;
        std::cout << "  Gradient + GB energy time: " << printMs(energyProfiler.grainGradientNs.load()) << " ms" << std::endl;
        std::cout << "  GB energy call time: " << printMs(energyProfiler.grainGbEnergyNs.load()) << " ms" << std::endl;
        std::cout << "  Grain energy time: " << printMs(energyProfiler.grainEnergyNs.load()) << " ms" << std::endl;
        std::cout << "  Total grain diff time: " << printMs(energyProfiler.grainTotalNs.load()) << " ms" << std::endl;
    }
    if (energyProfiler.particleCalls > 0) {
        std::cout << "Particle comp diff calls: " << energyProfiler.particleCalls << std::endl;
        std::cout << "  μ calculation time: " << printMs(energyProfiler.particleMuNs.load()) << " ms" << std::endl;
        std::cout << "  Total particle diff time: " << printMs(energyProfiler.particleTotalNs.load()) << " ms" << std::endl;
    }
    std::cout << "=== End Energy Profiling Results ===\n" << std::endl;
}

double calcGrainBoundaryEnergy(eulerAngles orient, const std::array<double,3>& grad3) {
    // Parameters
    const double MIN_GB = 2.9297;
    const double MAX_GB = 8.0;
    const double PI = 3.141592653589793;

    std::array<double,3> gradNormal = { -grad3[0], -grad3[1], 1.0 };

    double r = std::sqrt(gradNormal[0]*gradNormal[0] + gradNormal[1]*gradNormal[1] + gradNormal[2]*gradNormal[2]);
    if (r < 1e-12) return MIN_GB;

    gradNormal[0] /= r; gradNormal[1] /= r; gradNormal[2] /= r;

    std::array<double,3> surfPlaneVec = eulerRotate(orient, gradNormal);

    double x = surfPlaneVec[0], y = surfPlaneVec[1], z = surfPlaneVec[2];
    double mag = std::sqrt(x*x + y*y + z*z);
    if (mag < 1e-12) return MIN_GB;
    
    double zOverR = z / mag;
    if (zOverR > 1.0) zOverR = 1.0;
    if (zOverR < -1.0) zOverR = -1.0;
    double theta = std::acos(zOverR);
    double phi = std::atan2(y, x);
    if (!std::isfinite(theta) || !std::isfinite(phi)) return MIN_GB;

    if (SH_energyLUT.empty()) {
        std::cerr << "SH energy lookup table is empty. Please call readCoeffs() before running the simulation." << std::endl;
        return MIN_GB;
    }

    if (phi < 0.0) phi += 2.0 * PI;
    int thetaIndex = static_cast<int>((theta / PI) * (SH_LUT_THETA_SIZE - 1) + 0.5);
    int phiIndex = static_cast<int>((phi / (2.0 * PI)) * SH_LUT_PHI_SIZE);

    if (thetaIndex < 0) thetaIndex = 0;
    if (thetaIndex >= SH_LUT_THETA_SIZE) thetaIndex = SH_LUT_THETA_SIZE - 1;
    if (phiIndex < 0) phiIndex = 0;
    if (phiIndex >= SH_LUT_PHI_SIZE) phiIndex = SH_LUT_PHI_SIZE - 1;

    float energy = SH_energyLUT[phiIndex + thetaIndex * SH_LUT_PHI_SIZE];
    return static_cast<double>(energy);
}

std::array<double,3> eulerRotate(eulerAngles orient, std::array<double,3> rotatedVector) {
    double c1 = std::cos(orient.theta1), s1 = std::sin(orient.theta1);
    double cP = std::cos(orient.phi),  sP = std::sin(orient.phi);
    double c2 = std::cos(orient.theta2), s2 = std::sin(orient.theta2);

    double R00 =  c1*c2 - s1*s2*cP;
    double R01 = -c1*s2 - s1*c2*cP;
    double R02 =  s1*sP;

    double R10 =  s1*c2 + c1*s2*cP;
    double R11 = -s1*s2 + c1*c2*cP;
    double R12 = -c1*sP;

    double R20 =  s2*sP;
    double R21 =  c2*sP;
    double R22 =  cP;

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

    double denom = mag1*mag2;
    double cosval = 0.0;
    if (denom > 1e-15) cosval = adotb / denom;
    else cosval = 1.0;

    if (cosval > 1.0) cosval = 1.0;
    if (cosval < -1.0) cosval = -1.0;

    double angle = acos(cosval);

    if (angle > 3.141592/2) angle = 3.141592/2;
    return angle;
}

float calcPhaseDiffEnergy(int nodeIdx, config mConfig) {
    auto& nd = globalField.nodes;
    energyProfiler.phaseCalls++;
    auto totalStart = std::chrono::steady_clock::now();

    auto thermoStart = std::chrono::steady_clock::now();
    double pha = nd.phase[nodeIdx];
    double part = nd.particleComp[nodeIdx];
    double uc = underCool(nd.temp[nodeIdx], mConfig.meltTemp);
    double driveToLiq = (0.0349*nd.temp[nodeIdx] - 101.0704)/mConfig.molarVolume;
    double driveToSol = (-0.0212 *nd.temp[nodeIdx] + 61.3952)/mConfig.molarVolume;
    double eLoc = (-2 + 2*pha) * uc * driveToSol  + (2 * pha * (1 - uc)*driveToLiq );
    eLoc = eLoc + mConfig.phaseCoefficient*(((2* pha-2) * pow(nd.maxGrainPhase[nodeIdx],2)));
    double sizeScale = 0.66*4*3.14149 * pow(mConfig.particleRadius,2)/((4/3)*3.14159*pow(mConfig.particleRadius,3));
    eLoc = eLoc +  (pow(part,2)*mConfig.particleSolidIntEnergy*sizeScale) - (pow(part,2)*mConfig.particleLiquidIntEnergy*sizeScale);
    auto thermoEnd = std::chrono::steady_clock::now();
    energyProfiler.phaseThermoNs += std::chrono::duration_cast<std::chrono::nanoseconds>(thermoEnd - thermoStart).count();

    auto lapStart = std::chrono::steady_clock::now();
    double sumPh = 0.0;
    double nCount = 0;
    for (int nb = 0; nb < 4; nb++) {
        int nbrIdx = nd.neighbors[nodeIdx].idx[nb];
        if (nbrIdx >= TOTAL_NODES) {
            sumPh += pha; // boundary condition
        } else {
            sumPh += nd.phase[nbrIdx];
        }
        nCount++;
    }
    double lapPh = 0.0;
    if (nCount > 0) {
        lapPh = sumPh - nCount * pha;
    }
    lapPh = -1*lapPh* (1/(mConfig.dx*mConfig.dx));
    double eGrad = lapPh * 0.5*mConfig.phaseGradCo;
    auto lapEnd = std::chrono::steady_clock::now();
    energyProfiler.phaseLaplacianNs += std::chrono::duration_cast<std::chrono::nanoseconds>(lapEnd - lapStart).count();

    double diffFree = eLoc + eGrad;
    diffFree = safeClamp(diffFree, -1e22, 1e22);
    auto totalEnd = std::chrono::steady_clock::now();
    energyProfiler.phaseTotalNs += std::chrono::duration_cast<std::chrono::nanoseconds>(totalEnd - totalStart).count();
    return (float)diffFree;
}

std::array<float,9> calcGrainDiffEnergy(int nodeIdx, config mConfig) {
    auto& nd = globalField.nodes;
    std::array<float,9> diffFree;
    energyProfiler.grainCalls++;
    auto totalStart = std::chrono::steady_clock::now();
    double sumOtherGrainsSquared = nd.sumGrains[nodeIdx];

    for (int g = 0; g < 9; g++) {
        double gra = nd.grainPhases[nodeIdx][g];
        int activeGrain = nd.activeGrains[nodeIdx][g];
        double pha = nd.phase[nodeIdx];
        if(activeGrain < 0) {
            diffFree[g] = 0.0;
            continue;
        }

        auto neighborStart = std::chrono::steady_clock::now();
        double grainLeft = 0.0, grainRight = 0.0, grainUp = 0.0, grainDown = 0.0;
        double grainUpLeft = 0.0, grainUpRight = 0.0, grainDownLeft = 0.0, grainDownRight = 0.0;
        for(int gg = 0; gg < 9; gg++) {
            int nbrIdx;
            nbrIdx = nd.neighbors[nodeIdx].idx[0];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainLeft = nd.grainPhases[nbrIdx][gg];
            nbrIdx = nd.neighbors[nodeIdx].idx[1];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainRight = nd.grainPhases[nbrIdx][gg];
            nbrIdx = nd.neighbors[nodeIdx].idx[2];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainUp = nd.grainPhases[nbrIdx][gg];
            nbrIdx = nd.neighbors[nodeIdx].idx[3];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainDown = nd.grainPhases[nbrIdx][gg];
            nbrIdx = nd.neighbors[nodeIdx].idx[4];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainUpLeft = nd.grainPhases[nbrIdx][gg];
            nbrIdx = nd.neighbors[nodeIdx].idx[5];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainUpRight = nd.grainPhases[nbrIdx][gg];
            nbrIdx = nd.neighbors[nodeIdx].idx[6];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainDownLeft = nd.grainPhases[nbrIdx][gg];
            nbrIdx = nd.neighbors[nodeIdx].idx[7];
            if(nbrIdx < TOTAL_NODES && activeGrain == nd.activeGrains[nbrIdx][gg])
                grainDownRight = nd.grainPhases[nbrIdx][gg];
        }
        auto neighborEnd = std::chrono::steady_clock::now();
        energyProfiler.grainNeighborNs += std::chrono::duration_cast<std::chrono::nanoseconds>(neighborEnd - neighborStart).count();

        auto countStart = std::chrono::steady_clock::now();
        double sumGr = grainLeft + grainRight + grainUp + grainDown + 0.7071*(grainUpLeft + grainUpRight + grainDownLeft + grainDownRight);
        double nCount = 0;
        for(int nb = 0; nb < 4; nb++) {
            int nbrIdx = nd.neighbors[nodeIdx].idx[nb];
            if(nbrIdx < TOTAL_NODES && nd.exists[nbrIdx] != 0) nCount++;
        }
        for(int nb = 4; nb < 8; nb++) {
            int nbrIdx = nd.neighbors[nodeIdx].idx[nb];
            if(nbrIdx < TOTAL_NODES && nd.exists[nbrIdx] != 0) nCount += 0.7071;
        }
        double lapGr = 0.0;
        if(nCount > 0) {
            lapGr = sumGr - nCount * gra;
        }
        double grainGrad = -1*lapGr * mConfig.grainGradCo / (mConfig.dx * mConfig.dx);
        auto countEnd = std::chrono::steady_clock::now();
        energyProfiler.grainCountNs += std::chrono::duration_cast<std::chrono::nanoseconds>(countEnd - countStart).count();

        auto gradientStart = std::chrono::steady_clock::now();
        double gx = 0.0, gy = 0.0, gz = 0.0;
        double inv2dx = 0.5 / mConfig.dx;
        gx = (grainRight - grainLeft) * inv2dx * 1.0 +
             0.5 * ( (grainUpRight - grainUpLeft) + (grainDownRight - grainDownLeft) ) * inv2dx;
        gy = (grainDown - grainUp) * inv2dx * 1.0 +
             0.5 * ( (grainDownRight - grainUpRight) + (grainDownLeft  - grainUpLeft) ) * inv2dx;
        std::array<double,3> localGrad3 = { gx, gy, gz };
        auto gbStart = std::chrono::steady_clock::now();
        double gbEnergy = calcGrainBoundaryEnergy(nd.orientations[nodeIdx][g], localGrad3);
        auto gbEnd = std::chrono::steady_clock::now();
        energyProfiler.grainGbEnergyNs += std::chrono::duration_cast<std::chrono::nanoseconds>(gbEnd - gbStart).count();
        if (gbEnergy < 2.92) {
            gbEnergy = 2.92;
        }
        grainGrad = 0.5 * grainGrad * gbEnergy;
        auto gradientEnd = std::chrono::steady_clock::now();
        energyProfiler.grainGradientNs += std::chrono::duration_cast<std::chrono::nanoseconds>(gradientEnd - gradientStart).count();

        auto energyStart = std::chrono::steady_clock::now();
        double comp = sumOtherGrainsSquared - (gra*gra);
        double grainEnergy = mConfig.grainPreCo*(15*(pow(gra,3)-gra)*underCool(nd.temp[nodeIdx], mConfig.meltTemp) + (500*gra*comp*mConfig.grainIntWidth));
        grainGrad = safeClamp(grainGrad, -1e12, 1e12);
        grainEnergy = safeClamp(grainEnergy, -1e12, 1e12);
        diffFree[g] = grainGrad + grainEnergy;
        auto energyEnd = std::chrono::steady_clock::now();
        energyProfiler.grainEnergyNs += std::chrono::duration_cast<std::chrono::nanoseconds>(energyEnd - energyStart).count();
    }

    auto totalEnd = std::chrono::steady_clock::now();
    energyProfiler.grainTotalNs += std::chrono::duration_cast<std::chrono::nanoseconds>(totalEnd - totalStart).count();
    return diffFree;
}

double calcTemp(int nodeIdx, config& modelConfig, int t) {
    auto& nd = globalField.nodes;
    return modelConfig.startTemp + ((nd.tempDistance[nodeIdx]*modelConfig.tGrad) - modelConfig.tGrad*modelConfig.coolingRate* (modelConfig.dt*t));
}

double calcParticleCompDiff(int nodeIdx, config& modelConfig) {
    auto& nd = globalField.nodes;
    energyProfiler.particleCalls++;
    auto totalStart = std::chrono::steady_clock::now();

    auto muStart = std::chrono::steady_clock::now();
    double c = nd.particleComp[nodeIdx];
    double pha = nd.phase[nodeIdx];
    const double Vp = (4.0/3.0) * 3.14159 * std::pow(modelConfig.particleRadius, 3);
    const double cellPackedVolume = modelConfig.dx*modelConfig.dx*modelConfig.thickness2D;
    double sizeScale = 0.6*4*3.14149 * pow(modelConfig.particleRadius,2)/((4/3)*3.14159*pow(modelConfig.particleRadius,3));
    double muLocal = 2 *c  *pha * pow(modelConfig.dx,2)*modelConfig.thickness2D * modelConfig.particleSolidIntEnergy*sizeScale
                   +2 *c * (1.0 - pha) * pow(modelConfig.dx,2)*modelConfig.thickness2D*modelConfig.particleLiquidIntEnergy*sizeScale
                   + 2e7*(pow(modelConfig.dx,2)*modelConfig.thickness2D*(std::sin(2*3.14159 * (c*cellPackedVolume) / Vp)));
    double muGrav = 9.81 * (modelConfig.particleDensity-modelConfig.density) *pow(modelConfig.dx,4) *nd.yPos[nodeIdx];
    if (c > 1) muLocal  = muLocal+ 5e-17*std::pow(10*(c-1),2);
    if (c < 0) muLocal = muLocal-5e-17*std::pow(10*(c),2);
    double mu = muLocal + muGrav;
    auto muEnd = std::chrono::steady_clock::now();
    energyProfiler.particleMuNs += std::chrono::duration_cast<std::chrono::nanoseconds>(muEnd - muStart).count();

    if (!std::isfinite(mu) || std::isnan(mu)) {
        std::cerr << "[WARN] calcParticleCompDiff: non-finite mu at node id=" << nd.id[nodeIdx]
                  << "  c=" << c << "  particleRadius=" << modelConfig.particleRadius << "\n";
        mu = 0.0;
    }

    mu = safeClamp(mu, -1e12, 1e12);
    auto totalEnd = std::chrono::steady_clock::now();
    energyProfiler.particleTotalNs += std::chrono::duration_cast<std::chrono::nanoseconds>(totalEnd - totalStart).count();
    return mu;
}

// Stub for compatibility - actual free energy calculation would need refactoring
std::vector<std::vector<double>> calcFreeEnergy(gridField& field, config modelConfig) {
    std::vector<std::vector<double>> result(TOTAL_NODES);
    return result;
}
