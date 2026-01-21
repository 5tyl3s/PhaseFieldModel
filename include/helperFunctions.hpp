#ifndef HELPER_FUNCTIONS_HPP
#define HELPER_FUNCTIONS_HPP

#include <string>
#include <vector>
#include "gridField.hpp"

// BMP writer - saves framebuffer data to BMP file
bool saveBMP(const std::string &path, int width, int height, const std::vector<unsigned char> &rgba);

// File I/O functions for saving simulation results
void saveGridData(gridField &globalField);

#endif // HELPER_FUNCTIONS_HPP
