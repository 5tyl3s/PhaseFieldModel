#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include "gridField.hpp"
#include "matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

class FieldVisualizer {
private:
    int updateInterval;
    int gridRows, gridCols;
    
public:
    FieldVisualizer(int rows, int cols, int interval = 100) 
        : gridRows(rows), gridCols(cols), updateInterval(interval) {}

    void visualizeFields(const gridField& field, int timestep) {
        // Only update every updateInterval steps
        if (timestep % updateInterval != 0 && timestep != 0) return;

        // Extract phase field
        std::vector<std::vector<double>> phaseField(gridRows, std::vector<double>(gridCols));
        for (int i = 0; i < gridRows; i++) {
            for (int j = 0; j < gridCols; j++) {
                phaseField[i][j] = field.grid[i][j].phase;
            }
        }

        // Extract grain field (dominant grain ID)
        std::vector<std::vector<double>> grainField(gridRows, std::vector<double>(gridCols));
        for (int i = 0; i < gridRows; i++) {
            for (int j = 0; j < gridCols; j++) {
                double maxGrain = -1.0;
                for (int g = 0; g < field.grid[i][j].grainsHere; g++) {
                    if (field.grid[i][j].grainPhases[g] > maxGrain) {
                        maxGrain = field.grid[i][j].grainPhases[g];
                    }
                }
                grainField[i][j] = (maxGrain >= 0.0) ? maxGrain : 0.0;
            }
        }

        // Extract temperature field
        std::vector<std::vector<double>> tempField(gridRows, std::vector<double>(gridCols));
        for (int i = 0; i < gridRows; i++) {
            for (int j = 0; j < gridCols; j++) {
                tempField[i][j] = field.grid[i][j].temp;
            }
        }

        // Create figure with 3 subplots
        plt::figure_size(1200, 400);
        
        // Phase field
        plt::subplot(1, 3, 1);
        plt::imshow(phaseField, {{"cmap", "viridis"}});
        plt::colorbar();
        plt::title("Phase Field (t=" + std::to_string(timestep) + ")");
        plt::xlabel("X");
        plt::ylabel("Y");

        // Grain field
        plt::subplot(1, 3, 2);
        plt::imshow(grainField, {{"cmap", "plasma"}});
        plt::colorbar();
        plt::title("Grain Field (Max Phase)");
        plt::xlabel("X");
        plt::ylabel("Y");

        // Temperature field
        plt::subplot(1, 3, 3);
        plt::imshow(tempField, {{"cmap", "hot"}});
        plt::colorbar();
        plt::title("Temperature Field");
        plt::xlabel("X");
        plt::ylabel("Y");

        plt::tight_layout();
        plt::pause(0.001); // Non-blocking pause to allow window interaction
    }

    void setUpdateInterval(int interval) {
        updateInterval = interval;
    }
};