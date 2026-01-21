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

auto safeClamp = [](double v, double low, double high){
    if (std::isnan(v) || std::isinf(v)) return (low+high)/2.0;
    return std::max(low, std::min(high, v));
};

bool readCoeffs(const std::string& filename) {
    std::ifstream file(filename);
    std::cout << "Reading Spherical Harmonic Coefficients from: " << filename << std::endl;
    if(!file.is_open()) return false;

    SH_coeffs.clear();
    std::cout << "Loading Coefficients..." << std::endl;
    double val;
    while(file >> val) SH_coeffs.push_back(val);
    std::cout << "Loaded " << SH_coeffs.size() << " Coefficients." << std::endl;
    file.close();
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

double sumOtherGrainsSquared(int notIndex,node Node, int numGrains) {
    double notThisGrain = 0;
    for (int gr = 0; gr < 9;gr++) {
        notThisGrain = notThisGrain + (Node.grainPhases[gr]*Node.grainPhases[gr]);
    }
    return notThisGrain-(Node.grainPhases[notIndex]*Node.grainPhases[notIndex]);
    
}






double underCool(double temp, double meltTemp) {
    //0 Above melt, 1 below
    //std::cout << "Temperature: " << temp << " UnderCooled:  " << 0.5*(1-(tanh(1000000*((temp/meltTemp)-1)))) << std::endl;
    return 0.5*(1-(tanh(1000000*((temp/meltTemp)-1))));


}

double calcGrainBoundaryEnergy(eulerAngles orient, const std::array<double,3>& grad3) {
    // Parameters
    const double MIN_GB = 2.9297;   // J/m^2 (your minimum)
    const double MAX_GB = 8.0;    // upper clamp to avoid blow-ups (tunable)
    const double PI = 3.141592653589793;

    // 1) Build surface normal from gradient
    // If gradient is (dphi/dx, dphi/dy, dphi/dz), a reasonable normal is [-gx, -gy, 1] or [-gx, -gy, 0] depending on convention.
    // Here we assume the surface normal is approximately (-gx, -gy, 1) to include out-of-plane direction.
    std::array<double,3> gradNormal = { -grad3[0], -grad3[1], -grad3[2] };

    // If gradient magnitude is tiny -> no well-defined normal (flat), return min
    double r = std::sqrt(gradNormal[0]*gradNormal[0] + gradNormal[1]*gradNormal[1] + gradNormal[2]*gradNormal[2]);
    if (r < 1e-12) return MIN_GB;

    // Normalize
    gradNormal[0] /= r; gradNormal[1] /= r; gradNormal[2] /= r;

    // 2) Rotate normal by grain orientation
    std::array<double,3> surfPlaneVec = eulerRotate(orient, gradNormal);

    // 3) Convert to spherical coordinates (theta = colatitude in [0, pi], phi in [-pi, pi])
    double x = surfPlaneVec[0], y = surfPlaneVec[1], z = surfPlaneVec[2];
    double mag = std::sqrt(x*x + y*y + z*z);
    if (mag < 1e-12) return MIN_GB; // defensive
    // clamp z/mag into [-1,1]
    double zOverR = z / mag;
    if (zOverR > 1.0) zOverR = 1.0;
    if (zOverR < -1.0) zOverR = -1.0;
    double theta = std::acos(zOverR);
    double phi = std::atan2(y, x);
    if (!std::isfinite(theta) || !std::isfinite(phi)) return MIN_GB;

    // 4) Validate SH coeff vector size: expected (Lmax+1)^2 coefficients if stored in m=-l..l ordering
    int expected_count = (Lmax + 1) * (Lmax + 1);
    if ((int)SH_coeffs.size() < expected_count) {
        // coefficient file inconsistent: return min
        std::cerr << "SH_coeffs size " << SH_coeffs.size() << " < expected " << expected_count << std::endl;
        return MIN_GB;
    }

    // 5) Evaluate SH using ordering m = -l .. +l (this is the common ordering and matches MATLAB's typical layout)
    double energy = 0.0;
    int idx = 0;
    for (int l = 0; l <= Lmax; ++l) {
        for (int m = -l; m <= l; ++m) {
            // guard
            if (idx >= (int)SH_coeffs.size()) { return MIN_GB; }
            energy += SH_coeffs[idx++] * Y_lm_real(l, m, theta, phi);
        }
    }

    if (!std::isfinite(energy)) return MIN_GB;

    // 6) Clamp to expected range to avoid nonphysical explosion
    if (energy < MIN_GB) energy = MIN_GB;
    if (energy > MAX_GB) energy = MAX_GB;

    return energy;
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
    

    double denom = mag1*mag2;
    double cosval = 0.0;
    if (denom > 1e-15) cosval = adotb / denom;
    else cosval = 1.0; // treat as parallel if magnitude is tiny

    // clamp:
    if (cosval > 1.0) cosval = 1.0;
    if (cosval < -1.0) cosval = -1.0;

    double angle = acos(cosval);

    // optional: normalize to desired range without weird while loops
    if (angle > 3.141592/2) angle = 3.141592/2;
    return angle;

}

double calcPhaseDiffEnergy(node* nd, config mConfig) {
    //std::cout << "Calculating Phase Diff Energy..." << std::endl;
    double pha = nd->phase;
    double part = nd->particleComp;
    double uc = underCool(nd->temp, mConfig.meltTemp);
    double driveToLiq = 1000*(0.0349*nd->temp - 101.0704)/mConfig.molarVolume; // J/m^3
    double driveToSol = 1000*(-0.0212 *nd->temp + 61.3952)/mConfig.molarVolume; // J/m^3
    //std::cout << "PhaseCoefficient is: " << mConfig.phaseCoefficient << std::endl;
    double eLoc = (
        (-2 + 2*pha) * uc * driveToSol  + (2 * pha * (1 - uc)*driveToLiq )
    );
    //std::cout << "Phase: " << pha <<  " uc: " << uc <<" Local Term1: " << eLoc << std::endl;
    
    //std::cout << grainSum << std::endl;

    eLoc = eLoc + mConfig.phaseCoefficient *10000*((2* pha-2) * nd->sumGrains)+(2*pha* (1 - nd->sumGrains));
    
   double sizeScale = 0.66*4*3.14149 * pow(mConfig.particleRadius,2)/((4/3)*3.14159*pow(mConfig.particleRadius,3)); // surface area / volume
     

    eLoc = eLoc +  ((pha) * pow(part,2)*mConfig.particleSolidIntEnergy*sizeScale) + ((2*pha-2) * (pow(part,2))*mConfig.particleLiquidIntEnergy*sizeScale);
    //std::cout << "After Grain and Particle Term: " << eLoc << std::endl;
    // Simple 4-point Laplacian (von Neumann): neighbors 0,1,2,3 = left,right,up,down
    double sumPh = 0.0;
    int nCount = 0;
    for (int nb = 0; nb < 4; nb++) {
        node* nbr = nd->neighbors[nb];
        if (!nbr) continue;
        if (nbr->exists == 0) continue;
        sumPh += nbr->phase;
        nCount++;
    }
    double lapPh = 0.0;
    if (nCount > 0) {
        lapPh = sumPh - nCount * pha;  // (sum of neighbors - 4*center)
    }
    lapPh *= (1/(mConfig.dx*mConfig.dx));
    double eGrad = lapPh * mConfig.phaseGradCo;
    //std::cout << "Phase Gradient Coefficient: " << mConfig.phaseGradCo << std::endl;
    double diffFree = eLoc + eGrad;
    //std::cout << "Phase Diff Energy: " << diffFree << " GradientTerm: " << eGrad << " Gradient:" << eGrad << std::endl ;
    diffFree = safeClamp(diffFree, -1e8, 1e8);
    return diffFree;
 
 }
 
 

 
std::array<double,9> calcGrainDiffEnergy(node* nd, config mConfig) {
    std::array<double,9> diffFree;
    std::array<double,8> grainGradVals {0,0,0,0,0,0,0,0};

    for (int g = 0; g < 9; g++) {
        double gra = nd->grainPhases[g];
        int activeGrain = nd->activeGrains[g];
        double pha = nd->phase;
        if(activeGrain < 0) {
            diffFree[g] = 0.0;
            continue; // skip inactive slot
        }

        // Gather neighbor grain phases for same active grain (4-point only)
        double grainLeft = 0.0, grainRight = 0.0, grainUp = 0.0, grainDown = 0.0, grainUpLeft = 0.0, grainUpRight = 0.0, grainDownLeft = 0.0, grainDownRight = 0.0;
        for(int gg = 0; gg < 9; gg++) {
            node* nbrL = nd->neighbors[0];
            if(nbrL && activeGrain == nbrL->activeGrains[gg]) grainLeft = nbrL->grainPhases[gg];
            node* nbrR = nd->neighbors[1];
            if(nbrR && activeGrain == nbrR->activeGrains[gg]) grainRight = nbrR->grainPhases[gg];
            node* nbrU = nd->neighbors[2];
            if(nbrU && activeGrain == nbrU->activeGrains[gg]) grainUp = nbrU->grainPhases[gg];
            node* nbrD = nd->neighbors[3];
            if(nbrD && activeGrain == nbrD->activeGrains[gg]) grainDown = nbrD->grainPhases[gg];
            node* nbrUL = nd->neighbors[4];
            if(nbrUL && activeGrain == nbrUL->activeGrains[gg]) grainUpLeft = nbrUL->grainPhases[gg];
            node* nbrUR = nd->neighbors[5];
            if(nbrUR && activeGrain == nbrUR->activeGrains[gg]) grainUpRight = nbrUR->grainPhases[gg];
            node* nbrDL = nd->neighbors[6];
            if(nbrDL && activeGrain == nbrDL->activeGrains[gg]) grainDownLeft = nbrDL->grainPhases[gg];
            node* nbrDR = nd->neighbors[7];
            if(nbrDR && activeGrain == nbrDR->activeGrains[gg]) grainDownRight = nbrDR->grainPhases[gg];
        }

        // Simple 4-point Laplacian (von Neumann)
        double sumGr = grainLeft + grainRight + grainUp + grainDown + 0.7071*(grainUpLeft + grainUpRight + grainDownLeft + grainDownRight);
        int nCount = 0;
        for(int nb = 0; nb < 4; nb++) {
            node* nbr = nd->neighbors[nb];
            if(nbr && nbr->exists != 0) nCount++;
        } 
        for(int nb = 4; nb < 8; nb++) {
            node* nbr = nd->neighbors[nb];
            if(nbr && nbr->exists != 0) nCount+= 0.5;
        } 
        double lapGr = 0.0;
        if(nCount > 0) {
            lapGr = sumGr - nCount * gra;  // (sum of neighbors - 4*center)
        }
        // Gradient term drives growth: negative Laplacian (center surrounded by grain) lowers energy
        // Multiply by -1 so that lapGr > 0 (neighbors > center) makes grainGrad negative (lowers energy)
        double grainGrad = -lapGr * mConfig.grainGradCo / (mConfig.dx * mConfig.dx);

       double gx = 0.0, gy = 0.0, gz = 0.0;
        double inv2dx = 0.5 / mConfig.dx;   // (right - left)/(2*dx)
        double inv2dy = 0.5 / mConfig.dx;   // using dx for both dimensions assuming square grid

        gx = (grainRight - grainLeft) * inv2dx * 1.0 +               // horizontal
        0.5 * ( (grainUpRight - grainUpLeft)                    // diagonals
            + (grainDownRight - grainDownLeft) ) * inv2dx;

        gy = (grainDown - grainUp) * inv2dy * 1.0 +                  // vertical
        0.5 * ( (grainDownRight - grainUpRight)                 // diagonals
            + (grainDownLeft  - grainUpLeft) ) * inv2dy;
        // If you have any out-of-plane neighbor or z-derivative, compute gz similarly; otherwise leave gz=0.
 
        // Pass full 3D gradient (physical units) to GB energy routine
        std::array<double,3> localGrad3 = { gx, gy, gz };
        double gbEnergy = calcGrainBoundaryEnergy(nd->orientations[g], localGrad3);

        // continue using gbEnergy as you already do
        if (gbEnergy < 2.92) {
            std::cout << "Warning: Calculated GB energy below minimum. Clamping to minimum value." << std::endl;
            gbEnergy = 2.92;
        }
        grainGrad = grainGrad * gbEnergy*3; // scale gradient by local GB energy
        //std::cout << "Grain Gradient: " << grainGrad << std::endl;
        // Compute interaction/comp terms
        double comp = sumOtherGrainsSquared(g, *nd, 9);
        double notUC = 1.0 - underCool(nd->temp, mConfig.meltTemp);
        double grainEnergy = mConfig.grainPreCo*((pow(gra,3)-gra)*underCool(nd->temp, mConfig.meltTemp) + (gra*pow(comp,2)*mConfig.grainIntWidth) );
        grainGrad = safeClamp(grainGrad, -1e12, 1e12);
        grainEnergy = safeClamp(grainEnergy, -1e12, 1e12);
        diffFree[g] = grainGrad + grainEnergy;
        //std::cout << " Diff Energy: " << diffFree[g] << " GrainEnergyTerm: " << grainEnergy << " GrainGradTerm: " << grainGrad << std::endl;
    }

    return diffFree;
}
 
 
 double calcTemp(node* nd, config& modelConfig, int t) {
     
     return modelConfig.startTemp + ((nd->heightPos)*modelConfig.tGrad*modelConfig.dx) - modelConfig.tGrad*modelConfig.coolingRate* (modelConfig.dt*t);
     
 }
 
 double calcParticleCompDiff(node* nd, config& modelConfig) {
     // Compute CHEMICAL POTENTIAL (μ) which drives particle motion
     // Particles flow down the chemical potential gradient: J = -M * grad(μ)
     // μ = ∂f/∂c where f is the free energy density
     
     double c = nd->particleComp;
     double pha = nd->phase;
     
     // Chemical potential from local energy density:
     // f_particle = 2*c*pha^2*E_solid + 2*c*(1-pha)^2*E_liquid
     // ∂f/∂c = 2*pha^2*E_solid + 2*(1-pha)^2*E_liquid
     double sizeScale = 0.66*4*3.14149 * pow(modelConfig.particleRadius,2)/((4/3)*3.14159*pow(modelConfig.particleRadius,3)); // surface area / volume
     double muLocal = 2.0 *c * pha * pha * modelConfig.particleSolidIntEnergy*sizeScale
                    + 2.0 *c * (1.0 - pha) * (1.0 - pha) * modelConfig.particleLiquidIntEnergy*sizeScale;
     
     // Add gradient (interfacial) contribution: particles also respond to composition gradients
     // This provides a smoothing penalty for sharp concentration jumps
     // Gradient energy = kappa * |grad(c)|^2, so ∂/∂c_i includes laplacian term
     const double kappaParticle = 1e-8;  // interface energy coefficient (tune as needed)
     double lapC = 0;
     double sumNeighC = 0.0;
 

     int nExist = 0;
     for (int nb = 0; nb < 4; nb++) {
         node* nbr = nd->neighbors[nb];
         if (!nbr) continue;
         if (nbr->exists == 0) continue;
         sumNeighC += nbr->particleComp;
         nExist++;
     }
     lapC = sumNeighC - nExist * c;  // (sum of neighbors - 4*center)
     lapC *= (1/(modelConfig.dx*modelConfig.dx));

     
    // Variational chemical potential: μ = ∂f/∂c - kappa * laplacian(c)
    double mu = muLocal;

    // Defensive checks: guard against NaN/Inf coming from bad inputs (e.g. uninitialized particle radius)
    if (!std::isfinite(mu) || std::isnan(mu)) {
        std::cerr << "[WARN] calcParticleCompDiff: non-finite mu at node id=" << nd->id
                 // << "  muLocal=" << muLocal << "  lapC=" << lapC << "  pha=" << pha
                  << "  c=" << c << "  particleRadius=" << modelConfig.particleRadius << "\n";
        for (int nb = 0; nb < 4; ++nb) {
            node* nbr = nd->neighbors[nb];
            if (!nbr) continue;
            std::cerr << "    nbr["<<nb<<"].id="<<(nbr->id) << " c=" << nbr->particleComp << " pha=" << nbr->phase << "\n";
        }
        // fallback to safe zero potential
        mu = 0.0;
    }

    // clamp to avoid extremely large driving forces
    mu = safeClamp(mu, -1e12, 1e12);
    return mu;
 }



