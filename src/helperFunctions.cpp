#include "helperFunctions.hpp"
#include <fstream>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <iomanip>

// Simple BMP writer (expects data from glReadPixels in GL_RGBA, GL_UNSIGNED_BYTE)
bool saveBMP(const std::string &path, int width, int height, const std::vector<unsigned char> &rgba) {
    if ((int)rgba.size() < width * height * 4) return false;
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open()) return false;

    // BMP headers
    int rowBytes = width * 3;
    int pad = (4 - (rowBytes % 4)) % 4;
    int dataSize = (rowBytes + pad) * height;
    uint32_t fileSize = 54 + dataSize;

    unsigned char fileHeader[14] = {
        'B','M',
        0,0,0,0, // file size
        0,0, // reserved
        0,0,
        54,0,0,0
    };
    fileHeader[2] = (unsigned char)(fileSize      & 0xFF);
    fileHeader[3] = (unsigned char)((fileSize>>8) & 0xFF);
    fileHeader[4] = (unsigned char)((fileSize>>16)& 0xFF);
    fileHeader[5] = (unsigned char)((fileSize>>24)& 0xFF);

    unsigned char infoHeader[40] = {0};
    infoHeader[0] = 40;
    infoHeader[4] = (unsigned char)(width & 0xFF);
    infoHeader[5] = (unsigned char)((width>>8) & 0xFF);
    infoHeader[6] = (unsigned char)((width>>16) & 0xFF);
    infoHeader[7] = (unsigned char)((width>>24) & 0xFF);
    infoHeader[8] = (unsigned char)(height & 0xFF);
    infoHeader[9] = (unsigned char)((height>>8) & 0xFF);
    infoHeader[10] = (unsigned char)((height>>16) & 0xFF);
    infoHeader[11] = (unsigned char)((height>>24) & 0xFF);
    infoHeader[12] = 1; // planes
    infoHeader[14] = 24; // bits per pixel

    f.write(reinterpret_cast<char*>(fileHeader), sizeof(fileHeader));
    f.write(reinterpret_cast<char*>(infoHeader), sizeof(infoHeader));

    // glReadPixels returns pixels bottom-to-top by default (y=0 bottom). BMP expects bottom-up so write rows in order.
    // Convert RGBA -> BGR and write per-row with padding.
    std::vector<unsigned char> rowbuf(rowBytes + pad);
    for (int y = 0; y < height; ++y) {
        int base = y * width * 4;
        unsigned char *p = rowbuf.data();
        for (int x = 0; x < width; ++x) {
            unsigned char r = rgba[base + x*4 + 0];
            unsigned char g = rgba[base + x*4 + 1];
            unsigned char b = rgba[base + x*4 + 2];
            *p++ = b;
            *p++ = g;
            *p++ = r;
        }
        // padding bytes already zero-initialized if any
        for (int k = 0; k < pad; ++k) *p++ = 0;
        f.write(reinterpret_cast<char*>(rowbuf.data()), rowbuf.size());
    }

    f.close();
    return true;
}

// Save all grid data to CSV files
void saveGridData(gridField &globalField) {
    // Write Temperature Grid
    std::string fileNameTemp = "TempGrid.csv";
    std::ofstream output_file1(fileNameTemp);
    if (!output_file1.is_open()) {
        std::cerr << "Error: Could not open file " << fileNameTemp << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            output_file1 << globalField.nodes.temp[idx];
            if (j < GRID_COLS - 1) output_file1 << ",";
        }
        output_file1 << "\n";
    }
    output_file1.close();

    // Write Phase Grid
    std::string fileNamePhase = "PhaseGrid.csv";
    std::ofstream output_file2(fileNamePhase);
    if (!output_file2.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePhase << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            output_file2 << globalField.nodes.phase[idx];
            if (j < GRID_COLS - 1) output_file2 << ",";
        }
        output_file2 << "\n";
    }
    output_file2.close();

    // Write Grain Grid
    std::string fileNameGrain = "GrainGrid.csv";
    std::ofstream output_file3(fileNameGrain);
    if (!output_file3.is_open()) {
        std::cerr << "Error: Could not open file " << fileNameGrain << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            double tempOut = 0;
            int gNum = 0;
            for (int g = 0; g < 9; g++) {
                if (globalField.nodes.grainPhases[idx][g] > tempOut) {
                    tempOut = globalField.nodes.grainPhases[idx][g];
                    gNum = globalField.nodes.activeGrains[idx][g] + 1;
                }
            }
            output_file3 << gNum;
            if (j < GRID_COLS - 1) {
                output_file3 << ',';
            }
        }
        output_file3 << std::endl;
    }
    output_file3.close();

    // Write Particle Grid
    std::string fileNamePart = "ParticleGrid.csv";
    std::ofstream output_file4(fileNamePart);
    if (!output_file4.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePart << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            output_file4 << globalField.nodes.particleComp[idx];
            if (j < GRID_COLS - 1) {
                output_file4 << ',';
            }
        }
        output_file4 << std::endl;
    }
    output_file4.close();

    // Write Grain Orientations (Euler angles + quaternion)
    int nGrains = globalField.numGrains;
    if (nGrains > 0) {
        std::vector<bool> found(nGrains, false);
        std::vector<eulerAngles> eulers(nGrains);
        std::vector<int> counts(nGrains, 0);

        for (int i = 0; i < GRID_ROWS; ++i) {
            for (int j = 0; j < GRID_COLS; ++j) {
                int idx = i * GRID_COLS + j;
                for (int s = 0; s < globalField.nodes.grainsHere[idx]; ++s) {
                    int gid = globalField.nodes.activeGrains[idx][s];
                    if (gid < 0 || gid >= nGrains) continue;
                    counts[gid]++;
                    if (!found[gid]) {
                        found[gid] = true;
                        eulers[gid] = globalField.nodes.orientations[idx][s];
                    }
                }
            }
        }

        std::ofstream outG("GrainOrientations.csv");
        if (!outG.is_open()) {
            std::cerr << "Error: Could not open GrainOrientations.csv for writing" << std::endl;
        } else {
            outG << "grain_id,theta1_deg,Phi_deg,theta2_deg,q0,q1,q2,q3,node_count\n";
            for (int gid = 0; gid < nGrains; ++gid) {
                if (!found[gid]) {
                    outG << gid << ",NaN,NaN,NaN,NaN,NaN,NaN,NaN,0\n";
                    continue;
                }
                double phi1 = eulers[gid].theta1 * M_PI / 180.0;
                double Phi = eulers[gid].phi * M_PI / 180.0;
                double phi2 = eulers[gid].theta2 * M_PI / 180.0;

                double c1 = cos(phi1), s1 = sin(phi1);
                double cP = cos(Phi), sP = sin(Phi);
                double c2 = cos(phi2), s2 = sin(phi2);

                double R00 = c1 * c2 - s1 * s2 * cP;
                double R01 = -c1 * s2 - s1 * c2 * cP;
                double R02 = s1 * sP;
                double R10 = s1 * c2 + c1 * s2 * cP;
                double R11 = -s1 * s2 + c1 * c2 * cP;
                double R12 = -c1 * sP;
                double R20 = s2 * sP;
                double R21 = c2 * sP;
                double R22 = cP;

                // Quaternion from rotation matrix (w, x, y, z)
                double tr = R00 + R11 + R22;
                double qw, qx, qy, qz;
                if (tr > 0.0) {
                    double S = std::sqrt(tr + 1.0) * 2.0; // S=4*qw
                    qw = 0.25 * S;
                    qx = (R21 - R12) / S;
                    qy = (R02 - R20) / S;
                    qz = (R10 - R01) / S;
                } else if ((R00 > R11) & (R00 > R22)) {
                    double S = std::sqrt(1.0 + R00 - R11 - R22) * 2.0; // S=4*qx
                    qw = (R21 - R12) / S;
                    qx = 0.25 * S;
                    qy = (R01 + R10) / S;
                    qz = (R02 + R20) / S;
                } else if (R11 > R22) {
                    double S = std::sqrt(1.0 + R11 - R00 - R22) * 2.0; // S=4*qy
                    qw = (R02 - R20) / S;
                    qx = (R01 + R10) / S;
                    qy = 0.25 * S;
                    qz = (R12 + R21) / S;
                } else {
                    double S = std::sqrt(1.0 + R22 - R00 - R11) * 2.0; // S=4*qz
                    qw = (R10 - R01) / S;
                    qx = (R02 + R20) / S;
                    qy = (R12 + R21) / S;
                    qz = 0.25 * S;
                }

                outG << gid << "," << eulers[gid].theta1 << "," << eulers[gid].phi << "," << eulers[gid].theta2 << ","
                     << qw << "," << qx << "," << qy << "," << qz << "," << counts[gid] << "\n";
            }
            outG.close();
        }
    }
    
    std::cout << "Calc Completed, Saved Data" << std::endl;
}
