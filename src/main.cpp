#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include "importConfig.hpp"
#include "gridField.hpp"
#include "glVisualization.hpp"   // <- add: declare GLVisualizer
#include <GL/glew.h>             // <- add: for glReadPixels / GL types
#include <GLFW/glfw3.h>          // <- add: for glfwGetFramebufferSize / context
#include <omp.h>
#include <direct.h> // For _getcwd on Windows

#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <filesystem>
#include <iomanip>

// Use compile-time TOTAL_NODES from header so sizes match the grid definition
std::array<double, TOTAL_NODES> newTemp;
std::array<std::array<double,9>, TOTAL_NODES> grainDiffEn;
std::array<double, TOTAL_NODES> phaseDiffEn;
std::array<double, TOTAL_NODES> tempGrad;
std::array<double, TOTAL_NODES> tempPartComp;

// Simple BMP writer (expects data from glReadPixels in GL_RGBA, GL_UNSIGNED_BYTE)
static bool saveBMP(const std::string &path, int width, int height, const std::vector<unsigned char> &rgba) {
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

int main() {
    std::cout << "Hello, World!" << std::endl;
    std::cout << "Attempting to read config file..." << std::endl << std::flush;
    if(!readCoeffs("../config/spherical_harmonics_coeffs.txt")) {
        std::cerr << "Failed to read spherical harmonic coefficients!\n";
        return 1;
    }
    config configData = inputConfig("../config/modelConfig.json");
    if (configData.success == 0) {
        std::cerr << "Config Failed to load from: ../config/modelConfig.json" << std::endl << std::flush;
        return 1;
    }

    //gridField model2;
    globalField.init(configData);
    char cwd2[256];
    if (_getcwd(cwd2, sizeof(cwd2)) != NULL) {
        std::cout << "Current working directory: " << cwd2 << std::endl << std::flush;
    } else {
        std::cerr << "Could not get current working directory" << std::endl << std::flush;
    }    std::cout << "About to start PhaseField Sim..." << std::endl << std::flush;
    auto start1 = std::chrono::steady_clock::now();
    node* nd2;
    std::cout << "Starting Time Steps..." << std::endl << std::flush;

    // create visualizer before the simulation so it updates live
    GLVisualizer visualizer(GRID_ROWS, GRID_COLS, 50);
    bool haveViz = false;
    if (visualizer.initialize("Phase Field Visualization")) {
        haveViz = true;
        visualizer.updateFromGrid(globalField);
    } else {
        std::cerr << "Visualizer init failed, continuing without live display\n";
    }

    // prepare frames output dir
    std::filesystem::path framesDir = std::filesystem::current_path() / "viz_frames";
    std::error_code ec;
    std::filesystem::create_directories(framesDir, ec);

    int frameCounter = 0;
    int saveEvery = visualizer.getUpdateInterval();
    bool simulationComplete = false;

    for (int t = 0; t < configData.timeSteps; t++) {
        if (t%10 == 0) std::cout << t << "/" << configData.timeSteps << std::endl;
        
        // Check every 1000 timesteps if fully solidified with all grains present
        if (t > 0 && (t % 1000 == 0)) {
            bool allSolidified = true;
            bool allHaveGrains = true;
            for (int i = 0; i < GRID_ROWS && allSolidified; ++i) {
                for (int j = 0; j < GRID_COLS && allSolidified; ++j) {
                    node& n = globalField.grid[i][j];
                    // Check if solidified (phase >= 1.0)
                    if (n.phase < 0.9999) {
                        allSolidified = false;
                    }
                    // Check if has a grain
                    if (n.grainsHere <= 0) {
                        allHaveGrains = false;
                    }
                }
            }
            if (allSolidified && allHaveGrains) {
                std::cout << "Fully solidified with grains at timestep " << t << ", ending simulation\n";
                simulationComplete = true;
                break;
            }
        }
        #pragma omp parallel for private(nd2)
        for (int node = 0; node < TOTAL_NODES; node++) {
            nd2 = globalField.allNodes[node];
            tempGrad[node] = calcTemp(nd2,configData,t);
            phaseDiffEn[node] = calcPhaseDiffEnergy(nd2, configData);
            tempPartComp[node] = calcParticleCompDiff(nd2, configData);
            grainDiffEn[node] = calcGrainDiffEnergy(nd2, configData);
        }
        

        globalField.update(phaseDiffEn, grainDiffEn, tempPartComp, tempGrad);

        // update visualizer and save a frame periodically
        if (haveViz && (t % saveEvery == 0)) {
            visualizer.updateFromGrid(globalField);
            visualizer.processInput();
            visualizer.render();

            // read back full framebuffer and save BMP
            int winW = 0, winH = 0;
            // get framebuffer size from current GLFW context
            GLFWwindow* ctx = glfwGetCurrentContext();
            if (ctx) {
                glfwGetFramebufferSize(ctx, &winW, &winH);
            } else {
                winW = 0; winH = 0;
            }
            // fallback: use GRID_COLS/ROWS scaled to pixels; assume 800x800 if 0
            if (winW == 0 || winH == 0) { winW = GRID_COLS; winH = GRID_ROWS; }

            // allocate buffer and read pixels
            std::vector<unsigned char> pixels(winW * winH * 4);
            glReadPixels(0, 0, winW, winH, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());

            std::ostringstream fname;
            fname << framesDir.string() << "/frame_" << std::setw(6) << std::setfill('0') << frameCounter << ".bmp";
            if (!saveBMP(fname.str(), winW, winH, pixels)) {
                std::cerr << "Failed to write frame " << fname.str() << std::endl;
            } else {
                frameCounter++;
            }
        }

        // Stop simulation if visualizer window is closed by user
        if (haveViz && visualizer.shouldClose()) {
            std::cout << "Visualizer window closed by user, stopping simulation\n";
            break;
        }
    }

    // If simulation completed naturally or by solidification check, pause visualizer on last frame
    if (simulationComplete && haveViz) {
        std::cout << "Simulation complete. Displaying final frame in visualizer window.\n";
        std::cout << "Close the window to finish and save output files.\n";
        visualizer.updateFromGrid(globalField);
        visualizer.render();
        // Keep displaying the final frame until user closes window
        while (haveViz && !visualizer.shouldClose()) {
            visualizer.processInput();
            visualizer.render();
            std::this_thread::sleep_for(std::chrono::milliseconds(16)); // ~60 FPS
        }
    }

    // Write output files using GRID_ROWS / GRID_COLS
    std::string fileNameTemp = "TempGrid.csv";
    std::ofstream output_file1(fileNameTemp);
    if (!output_file1.is_open()) {
        std::cerr << "Error: Could not open file " << fileNameTemp << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            output_file1 << /* value */ globalField.grid[i][j].temp;
            if (j < GRID_COLS - 1) output_file1 << ",";
        }
        output_file1 << "\n";
    }
    output_file1.close();
    std::string fileNamePhase = "PhaseGrid.csv";
    std::ofstream output_file2(fileNamePhase);
    if (!output_file2.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePhase << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            int idx = i * GRID_COLS + j;
            output_file2  << globalField.grid[i][j].phase;
            if (j < GRID_COLS - 1) output_file2  << ",";
        }
        output_file2  << "\n";
    }
    std::cout << "Calc Completed, Saved Data";
    output_file2.close();

    std::string fileNameGrain = "GrainGrid.csv";
    std::ofstream output_file3(fileNameGrain);
    double tempOut = 0;
    int gNum;
    
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) { 
            gNum = 0;
            tempOut = 0;
            for (int g = 0; g < 9; g++) {
                if (globalField.grid[i][j].grainPhases[g] > tempOut) {
                    tempOut = globalField.grid[i][j].grainPhases[g];
                    gNum = globalField.grid[i][j].activeGrains[g]+1;
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


    std::string fileNamePart = "ParticleGrid.csv";
    std::ofstream output_file4(fileNamePart);
    if (!output_file4.is_open()) {
        std::cerr << "Error: Could not open file " << fileNamePart << std::endl;
    }
    for (int i = 0; i < GRID_ROWS; i++) {
        for (int j = 0; j < GRID_COLS; j++) {
            output_file4  << globalField.grid[i][j].particleComp;
            if (j < GRID_COLS - 1) {
                output_file4  << ',';
            }
        }
        output_file4  << std::endl;
    }
    std::cout << "Calc Completed, Saved Data";
    output_file4   .close();

    // --- Write grain orientations (Euler angles + quaternion) ---
    int nGrains = globalField.numGrains;
    if (nGrains > 0) {
        std::vector<bool> found(nGrains, false);
        std::vector<eulerAngles> eulers(nGrains);
        std::vector<int> counts(nGrains, 0);

        for (int i = 0; i < GRID_ROWS; ++i) {
            for (int j = 0; j < GRID_COLS; ++j) {
                node &nd = globalField.grid[i][j];
                for (int s = 0; s < nd.grainsHere; ++s) {
                    int gid = nd.activeGrains[s];
                    if (gid < 0 || gid >= nGrains) continue;
                    counts[gid]++;
                    if (!found[gid]) {
                        found[gid] = true;
                        eulers[gid] = nd.orientations[s];
                    }
                }
            }
        }

        std::ofstream outG("GrainOrientations.csv");
        if (!outG.is_open()) std::cerr << "Error: Could not open GrainOrientations.csv for writing" << std::endl;
        else {
            outG << "grain_id,theta1_deg,Phi_deg,theta2_deg,q0,q1,q2,q3,node_count\n";
            for (int gid = 0; gid < nGrains; ++gid) {
                if (!found[gid]) {
                    outG << gid << ",NaN,NaN,NaN,NaN,NaN,NaN,NaN,0\n";
                    continue;
                }
                double phi1 = eulers[gid].theta1 * M_PI / 180.0;
                double Phi  = eulers[gid].phi    * M_PI / 180.0;
                double phi2 = eulers[gid].theta2 * M_PI / 180.0;

                double c1 = cos(phi1), s1 = sin(phi1);
                double cP = cos(Phi),  sP = sin(Phi);
                double c2 = cos(phi2), s2 = sin(phi2);

                double R00 = c1*c2 - s1*s2*cP;
                double R01 = -c1*s2 - s1*c2*cP;
                double R02 = s1*sP;
                double R10 = s1*c2 + c1*s2*cP;
                double R11 = -s1*s2 + c1*c2*cP;
                double R12 = -c1*sP;
                double R20 = s2*sP;
                double R21 = c2*sP;
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

    std::cout << "Elapsed(ms)=" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start1).count() << std::endl;
    return 0;
}

// Helper function to resize grainDiffEn everywhere

