#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <string>
#include <vector>
#include <array>
#include "gridField.hpp"

// GLVisualizer: simple OpenGL viewer for phase/grain/temp/particle fields
class GLVisualizer {
private:
    GLFWwindow* window;
    GLuint VAO, VBO, EBO;
    GLuint shaderProgram;
    GLuint texturePhase, textureGrain, textureTemp, textureParticle;
    
    int gridRows, gridCols;
    int currentTimestep;
    int totalTimesteps;
    int updateInterval;
    bool isPaused;
    int displayMode; // 0: phase, 1: grain, 2: temp, 3: particle

    bool initialized_ = false;
    
    // optional CSV-backed storage (not required for live update)
    std::vector<std::vector<std::vector<float>>> phaseData;
    std::vector<std::vector<std::vector<float>>> grainData;
    std::vector<std::vector<std::vector<float>>> tempData;
    std::vector<std::vector<std::vector<float>>> particleData;
    
public:
    GLVisualizer(int rows, int cols, int interval = 100);
    ~GLVisualizer();
    bool initialize(const std::string& title = "Phase Field Visualization");
    void loadCSVData(const std::string& phaseFile, const std::string& grainFile,
                     const std::string& tempFile, const std::string& particleFile);
    void updateFromGrid(const gridField& field);
    void render();
    bool shouldClose();
    void processInput();
    void setUpdateInterval(int interval);
    void togglePause();
    void nextFrame();
    void prevFrame();
    void setDisplayMode(int mode); // 0-3

    int getUpdateInterval() const { return updateInterval; }
    bool isInitialized() const { return initialized_; }
    
private:
    GLuint createShaderProgram();
    void setupGeometry();
    void updateTexture(GLuint textureID, const std::vector<std::vector<float>>& data);
    void updateTextureTemp(GLuint textureID, const std::vector<std::vector<float>>& data, float startTemp);
    void updateTextureBlueYellow(GLuint textureID, const std::vector<std::vector<float>>& data);
    void updateTextureGrain(GLuint textureID, const std::vector<std::vector<float>>& data);
    void loadCSV(const std::string& filepath, std::vector<std::vector<std::vector<float>>>& out);
    std::string getShaderInfoLog(GLuint shader);
    std::string getProgramInfoLog(GLuint program);
};