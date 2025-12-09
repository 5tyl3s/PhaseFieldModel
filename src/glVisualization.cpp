#include "glVisualization.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <random>

// helper: pastel color generator (deterministic per grain id)
static std::unordered_map<int, std::array<unsigned char,3>> grainColorMap;
static std::array<unsigned char,3> pastelForID(int id) {
    if (id < 0) return {0,0,0};
    auto it = grainColorMap.find(id);
    if (it != grainColorMap.end()) return it->second;
    std::mt19937 rng(static_cast<unsigned int>(id * 1664525u + 1013904223u));
    std::uniform_real_distribution<float> d(0.5f, 0.95f); // pastel range
    float r = d(rng);
    float g = d(rng);
    float b = d(rng);
    float mix = 0.6f;
    float avg = (r+g+b)/3.0f;
    r = r*mix + avg*(1.0f-mix);
    g = g*mix + avg*(1.0f-mix);
    b = b*mix + avg*(1.0f-mix);
    std::array<unsigned char,3> col = {
        static_cast<unsigned char>(std::round(r*255.0f)),
        static_cast<unsigned char>(std::round(g*255.0f)),
        static_cast<unsigned char>(std::round(b*255.0f))
    };
    grainColorMap[id] = col;
    return col;
}

// Vertex Shader
const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec2 texCoord;

out vec2 TexCoord;

uniform mat4 projection;
uniform mat4 view;

void main()
{
    gl_Position = projection * view * vec4(position, 1.0);
    TexCoord = texCoord;
}
)";

// Fragment Shader
const char* fragmentShaderSource = R"(
#version 330 core
in vec2 TexCoord;
out vec4 FragColor;

uniform sampler2D texture1;

void main()
{
    vec4 texColor = texture(texture1, TexCoord);
    FragColor = texColor;
}
)";

GLVisualizer::GLVisualizer(int rows, int cols, int interval)
    : window(nullptr), gridRows(rows), gridCols(cols), 
      currentTimestep(0), totalTimesteps(1), updateInterval(interval),
      isPaused(false), displayMode(0)
{
}

GLVisualizer::~GLVisualizer() {
    if (window) {
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
        glDeleteVertexArrays(1, &VAO);
        glDeleteTextures(1, &texturePhase);
        glDeleteTextures(1, &textureGrain);
        glDeleteTextures(1, &textureTemp);
        glDeleteTextures(1, &textureParticle);
        glDeleteProgram(shaderProgram);
        glfwDestroyWindow(window);
        glfwTerminate();
    }
}

std::string GLVisualizer::getShaderInfoLog(GLuint shader) {
    int infoLogLength = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
    std::string infoLog(infoLogLength, ' ');
    glGetShaderInfoLog(shader, infoLogLength, nullptr, &infoLog[0]);
    return infoLog;
}

std::string GLVisualizer::getProgramInfoLog(GLuint program) {
    int infoLogLength = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);
    std::string infoLog(infoLogLength, ' ');
    glGetProgramInfoLog(program, infoLogLength, nullptr, &infoLog[0]);
    return infoLog;
}

GLuint GLVisualizer::createShaderProgram() {
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);
    
    int success;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        std::cerr << "Vertex shader compilation failed:\n" << getShaderInfoLog(vertexShader) << std::endl;
    }

    // Compile fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);
    
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        std::cerr << "Fragment shader compilation failed:\n" << getShaderInfoLog(fragmentShader) << std::endl;
    }

    // Link program
    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        std::cerr << "Shader program linking failed:\n" << getProgramInfoLog(program) << std::endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    
    return program;
}

void GLVisualizer::setupGeometry() {
    // Quad vertices (normalized device coordinates)
    float vertices[] = {
        // positions          // texture coords
         1.0f,  1.0f, 0.0f,   1.0f, 1.0f,
         1.0f, -1.0f, 0.0f,   1.0f, 0.0f,
        -1.0f, -1.0f, 0.0f,   0.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,   0.0f, 1.0f
    };
    
    unsigned int indices[] = {
        0, 1, 3,
        1, 2, 3
    };

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // texture coord attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

bool GLVisualizer::initialize(const std::string& title) {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return false;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(1400, 500, title.c_str(), nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return false;
    }

    glViewport(0, 0, 1400, 500);
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);

    shaderProgram = createShaderProgram();
    setupGeometry();

    // Create textures
    glGenTextures(1, &texturePhase);
    glGenTextures(1, &textureGrain);
    glGenTextures(1, &textureTemp);
    glGenTextures(1, &textureParticle);

    // ensure shader is active and set default uniforms (identity transforms)
    glUseProgram(shaderProgram);
    // set projection and view to identity so quad appears in NDC
    GLint projLoc = glGetUniformLocation(shaderProgram, "projection");
    GLint viewLoc = glGetUniformLocation(shaderProgram, "view");
    if (projLoc >= 0) {
        float identity[16] = {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        };
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, identity);
    }
    if (viewLoc >= 0) {
        float identity[16] = {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        };
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, identity);
    }
    // bind texture unit 0 to sampler once
    GLint texLoc = glGetUniformLocation(shaderProgram, "texture1");
    if (texLoc >= 0) glUniform1i(texLoc, 0);

    // mark initialized after successful setup
    initialized_ = true;
    return true;
}

void GLVisualizer::loadCSV(const std::string& filepath, std::vector<std::vector<std::vector<float>>>& out) {
    std::cout << "Attempting to load: " << filepath << std::endl;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open CSV file: " << filepath << std::endl;
        return;
    }

    std::string line;
    int lineCount = 0;
    while (std::getline(file, line)) {
        std::vector<std::vector<float>> frame;
        std::vector<float> row;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            row.push_back(std::stof(value));
        }

        std::cout << "Line " << lineCount << " has " << row.size() << " values" << std::endl;

        // Reshape row into 2D grid
        for (int i = 0; i < gridRows; i++) {
            std::vector<float> gridRow;
            for (int j = 0; j < gridCols; j++) {
                if (i * gridCols + j < row.size()) {
                    gridRow.push_back(row[i * gridCols + j]);
                } else {
                    gridRow.push_back(0.0f);
                }
            }
            frame.push_back(gridRow);
        }
        
        out.push_back(frame);
        lineCount++;
    }

    file.close();
    std::cout << "Loaded " << out.size() << " timesteps from " << filepath << std::endl;
}

void GLVisualizer::loadCSVData(const std::string& phaseFile, const std::string& grainFile,
                               const std::string& tempFile, const std::string& particleFile) {
    loadCSV(phaseFile, phaseData);
    loadCSV(grainFile, grainData);
    loadCSV(tempFile, tempData);
    loadCSV(particleFile, particleData);

    totalTimesteps = phaseData.size();
    std::cout << "Total timesteps: " << totalTimesteps << std::endl;
}

void GLVisualizer::updateTexture(GLuint textureID, const std::vector<std::vector<float>>& data) {
    // Convert normalized floats to 8-bit RGBA to avoid GL_FLOAT texture issues
    std::vector<unsigned char> pixelData;
    pixelData.reserve(gridRows * gridCols * 4);

    float minVal = 1e9f, maxVal = -1e9f;
    for (const auto& row : data) {
        for (float v : row) {
            if (v < minVal) minVal = v;
            if (v > maxVal) maxVal = v;
        }
    }
    if (maxVal == minVal) maxVal = minVal + 1.0f;
    float rangeInv = 1.0f / (maxVal - minVal);

    for (int i = 0; i < gridRows; ++i) {
        for (int j = 0; j < gridCols; ++j) {
            float normalized = (data[i][j] - minVal) * rangeInv;
            if (normalized < 0.0f) normalized = 0.0f;
            if (normalized > 1.0f) normalized = 1.0f;
            unsigned char c = static_cast<unsigned char>(normalized * 255.0f);
            pixelData.push_back(c);
            pixelData.push_back(c);
            pixelData.push_back(c);
            pixelData.push_back(255);
        }
    }

    glBindTexture(GL_TEXTURE_2D, textureID);
    // ensure byte alignment for tightly packed data
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    // upload as unsigned byte RGBA
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gridCols, gridRows, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelData.data());
    glBindTexture(GL_TEXTURE_2D, 0);
}

void GLVisualizer::updateTextureTemp(GLuint textureID, const std::vector<std::vector<float>>& data, float startTemp) {
    std::vector<unsigned char> pixelData;
    pixelData.reserve(gridRows * gridCols * 4);

    // Temperature scale: 300°C (573K) to 3500°C (3773K)
    const float minTempK = 2800.0f;   // 300°C
    const float maxTempK = 3000.0f;  // 3500°C
    const float meltTempK = 2896.0f; // melt temperature (white line)
    const float meltTolerance = 0.1f; // tolerance band around melt temp for white line

    for (int i = 0; i < gridRows; ++i) {
        for (int j = 0; j < gridCols; ++j) {
            float v = data[i][j];
            
            // Clamp to range
            if (v < minTempK) v = minTempK;
            if (v > maxTempK) v = maxTempK;
            
            // Normalize to [0, 1]
            float norm = (v - minTempK) / (maxTempK - minTempK);
            
            // Check if near melt temperature (draw white)
            if (std::abs(v - meltTempK) < meltTolerance) {
                pixelData.push_back(255);  // r
                pixelData.push_back(255);  // g
                pixelData.push_back(255);  // b
                pixelData.push_back(255);  // a
            } else {
                // Blue (cold, norm=0) -> Red (hot, norm=1)
                unsigned char r = static_cast<unsigned char>(std::round(255.0f * norm));
                unsigned char g = 0;
                unsigned char b = static_cast<unsigned char>(std::round(255.0f * (1.0f - norm)));
                pixelData.push_back(r);
                pixelData.push_back(g);
                pixelData.push_back(b);
                pixelData.push_back(255);
            }
        }
    }

    glBindTexture(GL_TEXTURE_2D, textureID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gridCols, gridRows, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelData.data());
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Phase & Particle: blue (0) -> yellow (1) mapping (r=v, g=v, b=1-v)
void GLVisualizer::updateTextureBlueYellow(GLuint textureID, const std::vector<std::vector<float>>& data) {
    std::vector<unsigned char> pixelData;
    pixelData.reserve(gridRows * gridCols * 4);

    float minVal = 0.0f;
    float maxVal = 1.0f;

    for (int i = 0; i < gridRows; ++i) {
        for (int j = 0; j < gridCols; ++j) {
            float v = data[i][j];
            float norm = (v - minVal) / (maxVal - minVal);
            if (norm < 0.0f) norm = 0.0f;
            if (norm > 1.0f) norm = 1.0f;
            unsigned char r = static_cast<unsigned char>(std::round(255.0f * norm));       // r = v
            unsigned char g = static_cast<unsigned char>(std::round(255.0f * norm));       // g = v
            unsigned char b = static_cast<unsigned char>(std::round(255.0f * (1.0f - norm))); // b = 1-v
            pixelData.push_back(r);
            pixelData.push_back(g);
            pixelData.push_back(b);
            pixelData.push_back(255);
        }
    }

    glBindTexture(GL_TEXTURE_2D, textureID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gridCols, gridRows, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelData.data());
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Grains: pastel random color per grain id (value contains id+1), 0 -> black / background
void GLVisualizer::updateTextureGrain(GLuint textureID, const std::vector<std::vector<float>>& data) {
    std::vector<unsigned char> pixelData;
    pixelData.reserve(gridRows * gridCols * 4);

    for (int i = 0; i < gridRows; ++i) {
        for (int j = 0; j < gridCols; ++j) {
            int val = static_cast<int>(std::round(data[i][j]));
            int id = val - 1; // stored as activeID+1 earlier
            if (id < 0) {
                pixelData.push_back(0);
                pixelData.push_back(0);
                pixelData.push_back(0);
                pixelData.push_back(255);
            } else {
                auto c = pastelForID(id);
                pixelData.push_back(c[0]);
                pixelData.push_back(c[1]);
                pixelData.push_back(c[2]);
                pixelData.push_back(255);
            }
        }
    }

    glBindTexture(GL_TEXTURE_2D, textureID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gridCols, gridRows, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelData.data());
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Replace calls in updateFromGrid to use specialized uploads
void GLVisualizer::updateFromGrid(const gridField& field) {
    // build temporary 2D float arrays for each field and update textures
    std::vector<std::vector<float>> phase(gridRows, std::vector<float>(gridCols));
    std::vector<std::vector<float>> temp(gridRows, std::vector<float>(gridCols));
    std::vector<std::vector<float>> particle(gridRows, std::vector<float>(gridCols));
    // grain texture will be built as RGBA bytes from orientations
    std::vector<unsigned char> grainPixels;
    grainPixels.reserve(gridRows * gridCols * 4);

    float maxTemp = -1e9f;

    for (int i = 0; i < gridRows; ++i) {
        for (int j = 0; j < gridCols; ++j) {
            const node& n = field.grid[i][j];
            phase[i][j] = static_cast<float>(n.phase);
            temp[i][j] = static_cast<float>(n.temp);
            particle[i][j] = static_cast<float>(n.particleComp);
            if (temp[i][j] > maxTemp) maxTemp = temp[i][j];
            // derive dominant grain slot and orientation
            float bestVal = -1.0f;
            int bestSlot = -1;
            for (int g = 0; g < n.grainsHere; ++g) {
                if (n.grainPhases[g] > bestVal) {
                    bestVal = static_cast<float>(n.grainPhases[g]);
                    bestSlot = g;
                }
            }
            // Build IPF color using cubic symmetry and the triangle defined by [100],[110],[111]
            if (bestSlot >= 0) {
                auto orient = n.orientations[bestSlot];
                // convert angles from degrees to radians (stored as degrees in gridField::addGrain)
                double phi1 = orient.theta1 * M_PI / 180.0;
                double Phi  = orient.phi    * M_PI / 180.0;
                double phi2 = orient.theta2 * M_PI / 180.0;

                // compute sample direction in crystal frame for sample Z (0,0,1)
                // use Bunge-like sequence: R = R_z(phi1) * R_x(Phi) * R_z(phi2)
                double c1 = std::cos(phi1), s1 = std::sin(phi1);
                double c2 = std::cos(phi2), s2 = std::sin(phi2);
                double cP = std::cos(Phi),  sP = std::sin(Phi);

                // full rotation matrix R (3x3)
                double R00 = c1*c2 - s1*cP*s2;
                double R01 = -c1*s2 - s1*cP*c2;
                double R02 = s1*sP;
                double R10 = s1*c2 + c1*cP*s2;
                double R11 = -s1*s2 + c1*cP*c2;
                double R12 = -c1*sP;
                double R20 = sP*s2;
                double R21 = sP*c2;
                double R22 = cP;

                // direction of sample z-axis in crystal coordinates (R^T * (0,0,1) == third column)
                double sx = R02;
                double sy = R12;
                double sz = R22;

                // Prepare cubic symmetry operators (24 proper rotations)
                std::vector<std::array<std::array<int,3>,3>> syms;
                int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
                for (int p = 0; p < 6; ++p) {
                    for (int s0 = -1; s0 <= 1; s0 += 2) {
                        for (int s1i = -1; s1i <= 1; s1i += 2) {
                            for (int s2 = -1; s2 <= 1; s2 += 2) {
                                std::array<std::array<int,3>,3> M = {{{0,0,0},{0,0,0},{0,0,0}}};
                                M[0][perms[p][0]] = s0;
                                M[1][perms[p][1]] = s1i;
                                M[2][perms[p][2]] = s2;
                                // compute determinant
                                int det =
                                    M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) -
                                    M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) +
                                    M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
                                if (det == 1) syms.push_back(M);
                            }
                        }
                    }
                }

                // choose symmetric variant that maps vector into fundamental triangle (maximize leading component)
                double bestA = -1.0;
                std::array<double,3> bestAbs = {0.0,0.0,0.0};
                for (const auto &M : syms) {
                    double tx = M[0][0]*sx + M[0][1]*sy + M[0][2]*sz;
                    double ty = M[1][0]*sx + M[1][1]*sy + M[1][2]*sz;
                    double tz = M[2][0]*sx + M[2][1]*sy + M[2][2]*sz;
                    double ax = std::abs(tx), ay = std::abs(ty), az = std::abs(tz);
                    // sort descending a>=b>=c
                    std::array<double,3> arr = {ax, ay, az};
                    std::sort(arr.begin(), arr.end(), std::greater<double>());
                    double a = arr[0];
                    if (a > bestA) {
                        bestA = a;
                        bestAbs = arr;
                    }
                }

                double a = bestAbs[0];
                double b = bestAbs[1];
                double c = bestAbs[2];
                if (a < 1e-12) {
                    // degenerate, paint black
                    grainPixels.push_back(0);
                    grainPixels.push_back(0);
                    grainPixels.push_back(0);
                    grainPixels.push_back(255);
                } else {
                    // barycentric-like weights for triangle vertices [100],[110],[111]
                    double w1 = (a - b);
                    double w2 = (b - c);
                    double w3 = c;
                    // normalize by a to get within triangle
                    w1 /= a; w2 /= a; w3 /= a;
                    // use weights as RGB and apply gamma for vibrancy
                    double gamma = 0.8;
                    double rc = std::pow(std::min(1.0, w1), gamma);
                    double gc = std::pow(std::min(1.0, w2), gamma);
                    double bc = std::pow(std::min(1.0, w3), gamma);
                    unsigned char rcb = static_cast<unsigned char>(std::round(rc * 255.0));
                    unsigned char gcb = static_cast<unsigned char>(std::round(gc * 255.0));
                    unsigned char bcb = static_cast<unsigned char>(std::round(bc * 255.0));
                    grainPixels.push_back(rcb);
                    grainPixels.push_back(gcb);
                    grainPixels.push_back(bcb);
                    grainPixels.push_back(255);
                }
            } else {
                // background
                grainPixels.push_back(0);
                grainPixels.push_back(0);
                grainPixels.push_back(0);
                grainPixels.push_back(255);
            }
        }
    }

    // Upload using requested colormaps
    updateTextureBlueYellow(texturePhase, phase);

    // Upload grain RGBA buffer generated from orientations
    glBindTexture(GL_TEXTURE_2D, textureGrain);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gridCols, gridRows, 0, GL_RGBA, GL_UNSIGNED_BYTE, grainPixels.data());
    glBindTexture(GL_TEXTURE_2D, 0);

    updateTextureTemp(textureTemp, temp, maxTemp);
    updateTextureBlueYellow(textureParticle, particle);
}

void GLVisualizer::render() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    int winW = 0, winH = 0;
    glfwGetFramebufferSize(window, &winW, &winH);
    if (winW == 0 || winH == 0) return;

    // cell sizes (no spacing) - exactly half the window
    int cellW = winW / 2;
    int cellH = winH / 2;

    std::array<GLuint,4> textures = { texturePhase, textureGrain, textureTemp, textureParticle };
    const char* labels[4] = { "Phase", "Grain", "Temperature", "Particle" };

    glUseProgram(shaderProgram);
    glActiveTexture(GL_TEXTURE0);

    // Draw each quad into its cell (no gaps)
    for (int idx = 0; idx < 4; ++idx) {
        int col = (idx % 2);      // 0 left, 1 right
        int row = (idx / 2);      // 0 top, 1 bottom

        int x = col * cellW;
        // OpenGL viewport origin is lower-left; row==0 is top => y = cellH
        int y = (row == 0) ? cellH : 0;

        glViewport(x, y, cellW, cellH);
        glBindTexture(GL_TEXTURE_2D, textures[idx]);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    }

    // restore full-window viewport
    glViewport(0, 0, winW, winH);

    // simple console status
    std::cout << "\rTimestep: " << currentTimestep << "/" << totalTimesteps
              << " | " << (isPaused ? "PAUSED" : "PLAYING") << "    ";
    std::cout.flush();

    glfwSwapBuffers(window);
}

bool GLVisualizer::shouldClose() {
    return glfwWindowShouldClose(window);
}

void GLVisualizer::processInput() {
    static bool spacePressed = false;
    
    glfwPollEvents();

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, true);
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
        if (!spacePressed) {
            togglePause();
            spacePressed = true;
        }
    } else {
        spacePressed = false;
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
        nextFrame();
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
        prevFrame();
    }
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) {
        setDisplayMode(0); // Phase
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) {
        setDisplayMode(1); // Grain
    }
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
        setDisplayMode(2); // Temperature
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS) {
        setDisplayMode(3); // Particle
    }

    // Auto-advance if not paused
    if (!isPaused) {
        currentTimestep++;
        if (currentTimestep >= totalTimesteps) {
            currentTimestep = 0;
        }
    }
}

void GLVisualizer::togglePause() {
    isPaused = !isPaused;
}

void GLVisualizer::nextFrame() {
    currentTimestep = (currentTimestep + 1) % totalTimesteps;
}

void GLVisualizer::prevFrame() {
    currentTimestep = (currentTimestep - 1 + totalTimesteps) % totalTimesteps;
}

void GLVisualizer::setDisplayMode(int mode) {
    if (mode >= 0 && mode <= 3) {
        displayMode = mode;
    }
}

void GLVisualizer::setUpdateInterval(int interval) {
    updateInterval = interval;
}