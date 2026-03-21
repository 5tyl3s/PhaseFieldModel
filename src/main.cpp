#include <iostream>
#include "model.hpp"
#include "modelStartup.hpp"
#include "importConfig.hpp"
#include "gridField.hpp"
#include "glVisualization.hpp"   // <- add: declare GLVisualizer
#include "helperFunctions.hpp"    // helper functions like saveBMP
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

// Use vectors for dynamic allocation (size determined at runtime)
std::vector<float> newTemp;
std::vector<std::array<float,9>> grainDiffEn;
std::vector<float> phaseDiffEn;
std::vector<float> tempGrad;
std::vector<float> tempPartComp;

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

    // Initialize grid with dimensions from config (default 500x500 if not specified)
    int gridRows = configData.gridRows;
    int gridCols = configData.gridCols;
    
    //gridField model2;
    globalField.init(configData, gridRows, gridCols);
    
    // Resize work vectors to match grid size
    newTemp.resize(globalField.totalNodes);
    grainDiffEn.resize(globalField.totalNodes);
    phaseDiffEn.resize(globalField.totalNodes);
    tempGrad.resize(globalField.totalNodes);
    tempPartComp.resize(globalField.totalNodes);
    
    auto cwd = std::filesystem::current_path();
    std::cout << "Current working directory: " << cwd << std::endl << std::flush;
    std::cout << "About to start PhaseField Sim..." << std::endl << std::flush;
    auto start1 = std::chrono::steady_clock::now();
    std::cout << "Starting Time Steps..." << std::endl << std::flush;

    // create visualizer before the simulation so it updates live
    GLVisualizer visualizer(globalField.gridRows, globalField.gridCols, 150);
    bool haveViz = false;
    if (configData.enableVisualization) {
        if (visualizer.initialize("Phase Field Visualization")) {
            haveViz = true;
            visualizer.updateFromGrid(globalField);
        } else {
            std::cerr << "Visualizer init failed, continuing without live display\n";
        }
    } else {
        std::cout << "Visualization disabled in config.\n";
    }

    // prepare frames output dir
    std::filesystem::path framesDir = std::filesystem::current_path() / "viz_frames";
    std::error_code ec;
    std::filesystem::create_directories(framesDir, ec);

    int frameCounter = 0;
    int saveEvery = visualizer.getUpdateInterval();
    bool simulationComplete = false;

    // Performance profiling
    long long totalStepTime = 0;
    long long totalVizTime = 0;
    long long totalFrameSaveTime = 0;
    long long totalSolidCheckTime = 0;
    long long totalOverheadTime = 0;
    int stepCount = 0;

    for (int t = 0; t < configData.timeSteps; t++) {
        auto stepStart = std::chrono::steady_clock::now();
        
        if (t%50 == 0) std::cout << t << "/" << configData.timeSteps << std::endl;
        
        // Check every 1000 timesteps if fully solidified with all grains present
        auto solidCheckStart = std::chrono::steady_clock::now();
        if (t > 0 && (t % 2000 == 0)) {
            bool allSolidified = true;
            bool allHaveGrains = true;
            for (int idx = 0; idx < globalField.totalNodes && allSolidified; ++idx) {
                // Check if solidified (phase >= 1.0)
                if (globalField.nodes.phase[idx] < 0.9999) {
                    allSolidified = false;
                }
                // Check if has a grain
                if (globalField.nodes.grainsHere[idx] <= 0) {
                    allHaveGrains = false;
                }
            }
            if (allSolidified && allHaveGrains) {
                std::cout << "Fully solidified with grains at timestep " << t << ", ending simulation\n";
                simulationComplete = true;
                break;
            }
        }
        auto solidCheckEnd = std::chrono::steady_clock::now();
        if (configData.enableProfiling) {
            totalSolidCheckTime += std::chrono::duration_cast<std::chrono::microseconds>(solidCheckEnd - solidCheckStart).count();
        }
        #pragma omp parallel for
        for (int node = 0; node < globalField.totalNodes; node++) {
            tempGrad[node] = calcTemp(node, configData, t);
            phaseDiffEn[node] = calcPhaseDiffEnergy(node, configData);
            tempPartComp[node] = calcParticleCompDiff(node, configData);
            grainDiffEn[node] = calcGrainDiffEnergy(node, configData);
        }
        auto diffEnergyEnd = std::chrono::steady_clock::now();
        long long diffEnergyTime = 0;
        if (configData.enableProfiling) {
            diffEnergyTime = std::chrono::duration_cast<std::chrono::microseconds>(diffEnergyEnd - stepStart).count();
        }

        globalField.update(phaseDiffEn, grainDiffEn, tempPartComp, tempGrad, configData.enableProfiling);
        if (configData.enableProfiling) {
            globalField.recordDiffEnergyTime(diffEnergyTime);
        }

        // update visualizer and save a frame periodically
        if (haveViz && (t % saveEvery == 0)) {
            auto vizStart = std::chrono::steady_clock::now();
            
            visualizer.updateFromGrid(globalField);
            visualizer.processInput();
            visualizer.render();

            auto frameSaveStart = std::chrono::steady_clock::now();
            
            // read back full framebuffer and save BMP
            int winW = 0, winH = 0;
            // get framebuffer size from current GLFW context
            GLFWwindow* ctx = glfwGetCurrentContext();
            if (ctx) {
                glfwGetFramebufferSize(ctx, &winW, &winH);
            } else {
                winW = 0; winH = 0;
            }
            // fallback: use gridCols/gridRows scaled to pixels; assume 800x800 if 0
            if (winW == 0 || winH == 0) { winW = globalField.gridCols; winH = globalField.gridRows; }

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
            
            auto frameSaveEnd = std::chrono::steady_clock::now();
            if (configData.enableProfiling) {
                totalVizTime += std::chrono::duration_cast<std::chrono::microseconds>(frameSaveStart - vizStart).count();
                totalFrameSaveTime += std::chrono::duration_cast<std::chrono::microseconds>(frameSaveEnd - frameSaveStart).count();
            }
        }

        // Stop simulation if visualizer window is closed by user
        if (haveViz && visualizer.shouldClose()) {
            std::cout << "Visualizer window closed by user, stopping simulation\n";
            break;
        }

        // Measure and average timestep performance
        if (configData.enableProfiling) {
            auto stepEnd = std::chrono::steady_clock::now();
            long long stepTime = std::chrono::duration_cast<std::chrono::microseconds>(stepEnd - stepStart).count();
            totalStepTime += stepTime;
            stepCount++;
        }

        if (configData.enableProfiling && stepCount >= 1000) {
            long long avgTime = totalStepTime / stepCount;
            long long avgVizTime = totalVizTime / stepCount;
            long long avgFrameSaveTime = totalFrameSaveTime / stepCount;
            long long avgSolidCheckTime = totalSolidCheckTime / stepCount;
            long long avgOverheadTime = avgTime - avgVizTime - avgFrameSaveTime - avgSolidCheckTime;
            
            std::cout << "\n=== Overall Timestep Performance (timesteps " << (t - 999) << "-" << t << ") ===" << std::endl;
            std::cout << "  Total Time per Step: " << avgTime << " μs" << std::endl;
            std::cout << "    - Calculation (from gridField): " << (avgTime - avgVizTime - avgFrameSaveTime - avgSolidCheckTime - avgOverheadTime) << " μs" << std::endl;
            std::cout << "    - Visualizer Update/Render: " << avgVizTime << " μs" << std::endl;
            std::cout << "    - Frame Save (glReadPixels+BMP): " << avgFrameSaveTime << " μs" << std::endl;
            std::cout << "    - Solidification Check: " << avgSolidCheckTime << " μs" << std::endl;
            std::cout << "    - Other Overhead: " << avgOverheadTime << " μs" << std::endl;
            std::cout << "======================================\n" << std::endl;
            
            totalStepTime = 0;
            totalVizTime = 0;
            totalFrameSaveTime = 0;
            totalSolidCheckTime = 0;
            totalOverheadTime = 0;
            stepCount = 0;
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

    // Save all grid data to files
    saveGridData(globalField);

    // Convert BMP frames to MP4 video at 30 FPS and delete BMP files
    std::cout << "\nConverting BMP frames to MP4 video...\n";
    std::string framesPattern = (framesDir / "frame_%06d.bmp").string();
    // Convert backslashes to forward slashes for ffmpeg compatibility
    std::replace(framesPattern.begin(), framesPattern.end(), '\\', '/');
    
    std::string outputVideo = (framesDir.parent_path() / "simulation_video.mp4").string();
    std::replace(outputVideo.begin(), outputVideo.end(), '\\', '/');
    
    // Build FFmpeg command using ffmpeg from PATH
    std::string ffmpegCmd = "ffmpeg -framerate 15 -i \"" + framesPattern + "\" -c:v mpeg4 -q:v 5 -y \"" + outputVideo + "\"";
    
    std::cout << "Running FFmpeg...\n";
    int ffmpegResult = system(ffmpegCmd.c_str());
    std::cout << "FFmpeg process exited with code: " << ffmpegResult << "\n";
    
    if (ffmpegResult == 0) {
        std::cout << "Video created successfully: " << outputVideo << "\n";
        
        // Delete all BMP files
        std::cout << "Deleting BMP frame files...\n";
        int deletedCount = 0;
        for (const auto& entry : std::filesystem::directory_iterator(framesDir)) {
            if (entry.path().extension() == ".bmp") {
                std::error_code ec;
                std::filesystem::remove(entry.path(), ec);
                if (!ec) {
                    deletedCount++;
                } else {
                    std::cerr << "Failed to delete: " << entry.path() << "\n";
                }
            }
        }
        std::cout << "Deleted " << deletedCount << " BMP files.\n";
    } else {
        std::cerr << "FFmpeg command failed with exit code: " << ffmpegResult << "\n";
        std::cerr << "BMP frames are still available in: " << framesDir << "\n";
        std::cerr << "Make sure all " << frameCounter << " frame BMP files are present.\n";
    }

    std::cout << "Elapsed(ms)=" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start1).count() << std::endl;
    return 0;
}


// Helper function to resize grainDiffEn everywhere

