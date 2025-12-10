#include <iostream>
#include "modelStartup.hpp"
#include <fstream>
#include <random>
#include <cmath>
#include <mutex>
#include <omp.h>
#include "gridField.hpp"

// modelStartup.cpp contains startup/config helpers and nucleation utilities.
// The main gridField::update implementation lives in src/gridField.cpp to
// keep the core update logic colocated with the grid data structures.

// If you need region-based or alternative update wrappers, add them here
// with distinct names (e.g. updateRegion(...)). Do NOT duplicate
// gridField::update(...) here to avoid multiple-definition/linker issues.








