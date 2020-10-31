/**
 * @file    statiff.cpp
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
using namespace statiff;
// using namespace cv;

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
    string  inputFileName    = ""; // command arg or config defined
    string outputFileName    = ""; // command arg or config defined
    // string outputFilePrefix  = ""; // none, output filenames will be the same as the standard
    // string outputFilePath    = ""; // same relative folder

    if (argInput) inputFileName = args::get(argInput); //input file is mandatory argument.
    if (inputFileName.empty()){ //not defined as command line argument? let's use config.yaml definition
        // ERROR! We do not have any definition of the input file
        logc.error ("main", "Input file missing. Please define it using --input='filename'");
        return ErrorCode::INVALID_ARG;
    }

    if (argInput) outputFileName = args::get(argOutput); //input file is mandatory argument.
    if (outputFileName.empty()){ //not defined as command line argument? let's use config.yaml definition
        // ERROR! We do not have any definition of the input file
        logc.error ("main", "Input file missing. Please define it using --output='filename'");
        return ErrorCode::INVALID_ARG;
    }

    int nThreads = DEFAULT_NTHREADS;
    if (argNThreads)   nThreads = args::get(argNThreads);

    // PRINT SUMMARY

    //LOAD GEOTIFF HEADER
    //EXTRACT PARAMETERS
    //CHECK PARAMETERS
    
    //PRINT SUMMARY

    //PRE-ALLOCATE VECTOR

    // FOR EACH ROW, COMPUTE:
        // Z <- READ PIXEL
        // CHECK IF DATA/NO_DATA
        // MIN, MAX, ACUM
        // B <-PARITTION BIN
        // UPDATE HIST(B)
    // NEXT

    // COMPUTE TOTAL AREA

    // FOR (Z:VECTOR)
        // DIST = Z - MEAN
        // VARIANCE_ACUM += DIST
    // NEXT

    // STDEV = SQRT(VARIANCE_ACUM)

    // SORT (V.BEING, V.END)
    // MEDIAN <- V.MIDDLE

    // DUMP TO OUTPUT FILE
    // PRINT SUMMARY ON SCREEN
    //------------>END

    return ErrorCode::NO_ERROR; // everything went smooth
}
