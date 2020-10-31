/**
 * @file    stat_tiff.cpp
 * @author  Jose Cappelletto (cappelletto@gmail.com)
 * @brief   geoTIFF statistics calculator tool. Use to analyze raster geoTIFF generated during the calculation of Landability and Measurability scores
 *          gdalinfo does not provide control over the histogram and the statistics calculation
 * @version 0.1
 * @date    2020-10-30
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "headers.h"
#include "helper.h"

#include "options.h"
#include "geotiff.hpp" // Geotiff class definitions

using namespace std;
using namespace cv;

logger::ConsoleOutput logc;

/*!
	@fn		int main(int argc, char* argv[])
	@brief	Main function
*/
int main(int argc, char *argv[])
{
    int retval = initParser(argc, argv);   // initial argument validation, populates arg parsing structure args
    if (retval != 0)  // some error ocurred, we have been signaled to stop
        return retval;
    std::ostringstream s;
    // Parameters hierarchy
    // ARGS > CONFIG > DEFAULT (this)
    // Input file priority: must be defined either by the config.yaml or --input argument
    string inputFileName    = ""; // command arg or config defined
    string outputFilePrefix = ""; // none, output filenames will be the same as the standard
    string outputFilePath   = ""; // same relative folder

    if (argInput) inputFileName = args::get(argInput); //input file is mandatory positional argument. Overrides any definition in configuration.yaml

    if (inputFileName.empty()){ //not defined as command line argument? let's use config.yaml definition
        // ERROR! We do not have any definition of the input file
        logc.error ("main", "Input file missing. Please define it using --input='filename' or inside a YAML configuration file (see --config option)");
        return -1;
    }

    // int nThreads = DEFAULT_NTHREADS;
    //     if (argNThreads)   nThreads = args::get(argNThreads);
    //     if (nThreads < 3) {
    //         s << "Number of used threads will be always 3 or higher. Asked for [" << yellow << nThreads << reset << "]" << endl;
    //         logc.warn("main", s);
    //     }
    // // override defaults or config file with command provided values (DEFAULT < CONFIG < ARGUMENT)
    // if (argAlphaRadius)     params.alphaShapeRadius = args::get(argAlphaRadius);
    // if (argGroundThreshold) params.groundThreshold  = args::get(argGroundThreshold);
    // if (argHeightThreshold) params.heightThreshold  = args::get(argHeightThreshold);
    // if (argSlopeThreshold)  params.slopeThreshold   = args::get(argSlopeThreshold);
    // if (argRobotHeight)     params.robotHeight      = args::get(argRobotHeight);
    // if (argRobotLength)     params.robotLength      = args::get(argRobotLength);
    // if (argRobotWidth)      params.robotWidth       = args::get(argRobotWidth);
    // if (argProtrusionSize)  params.protrusionSize   = args::get(argProtrusionSize);
    // if (argRotation){
    //                         params.rotation         = args::get(argRotation);
    //                         params.fixRotation      = true;
    // }   

}
