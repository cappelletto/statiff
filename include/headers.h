/**
 * @file headers.h
 * @author Jose Cappelletto (cappelletto@gmail.com)
 * @brief Single collection of global libraries required in major modules
 * @version 0.2
 * @date 2020-07-03
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef _PROJECT_HEADERS_H_

#define _PROJECT_HEADERS_H_
///Basic C and C++ libraries
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <omp.h>
/// OpenCV libraries. May need review for the final release
/// CUDA specific libraries
// #if USE_GPU
//     #include <opencv2/cudafilters.hpp>
//     #include "opencv2/cudafeatures2d.hpp"
//     #include "opencv2/xfeatures2d/cuda.hpp"
// #endif
#include "geotiff.hpp"
#include "statiff.hpp"

const std::string red("\033[1;31m");
const std::string green("\033[1;32m");
const std::string yellow("\033[1;33m");
const std::string blue("\033[1;34m");
const std::string purple("\033[1;35m");
const std::string cyan("\033[1;36m");

const std::string light_red("\033[0;31m");
const std::string light_green("\033[0;32m");
const std::string light_yellow("\033[0;33m");
const std::string light_blue("\033[0;34m");
const std::string light_purple("\033[0;35m");
const std::string light_cyan("\033[0;36m");

const std::string reset("\033[0m");
const std::string highlight("\033[30;43m");

#define DEFAULT_NTHREADS 12

#endif // _PROJECT_HEADERS_H_