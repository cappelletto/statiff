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
    double histMin = 0.0, histMax = 1.0;
    unsigned int nBins = 100;

    if (argNThreads) nThreads = args::get(argNThreads);
    if (argHistMax)   histMax = args::get(argHistMax);
    if (argHistMin)   histMin = args::get(argHistMin);
    if (argHistBin)     nBins = args::get(argHistBin);
    // validating input:
    if (histMin >= histMax){
        s << "Wrong histogram range arguments. histMin [" << yellow << histMin << reset << "] must be lower than histMax [" << yellow << histMax << reset << "]";
        logc.error("main", s);
        return statiff::INVALID_ARG;
    }
    
    if (nBins == 0){
        s << "Number of bins nBins[" << yellow << nBins << reset << "] must be positive integer";
        logc.error("main", s);
        return statiff::INVALID_ARG;
    }

    // PRINT SUMMARY
    logc.debug("main", "Summary information *****************************************************");
    cout << "\tInput file:   \t" << green << inputFileName << reset << endl;
    cout << "\tOutput file:  \t" << green << outputFileName << reset << endl;
    cout << "\tExport header:\t" << yellow << (argNoHeader ? "true" :  "false") << reset << endl;
    cout << "\tHistogram:    \t" << green << "[" << histMin << " : " << histMax << "]" << reset << "\t #Bins: " << yellow << nBins << reset << endl;
    cout << "\tNumber of threads (OpenMP):  \t" << yellow << nThreads << reset << endl;

    //LOAD GEOTIFF HEADER
    Geotiff inputGeotiff (inputFileName.c_str());
    if (!inputGeotiff.isValid()){ // check if nothing wrong happened with the constructor
        s << "Error opening geoTIFF file [" << inputFileName << "] via Geotiff class";

        return statiff::INVALID_DATA;
    }

    //**************************************
    // EXTRACT PARAMETERS FROM GEOTIFF
    // Get/print summary information of the TIFF 
    GDALDataset *poDataset;
    poDataset = inputGeotiff.GetDataset(); //pull the pointer to the main GDAL dataset structure
    // inputGeotiff.ShowInformation(); // show detailed info if asked for

    // canvas dimension
    int xSize = poDataset->GetRasterXSize();
    int ySize = poDataset->GetRasterYSize();
    // pixel resolution
    double xResolution, yResolution;
	double adfGeoTransform[6];
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
	{
	    xResolution = fabs(adfGeoTransform[1]);
	    yResolution = fabs(adfGeoTransform[5]);
	}
    // raster NO_DATA value
    int bGotNodata = 0;
    double dfNoData = GDALGetRasterNoDataValue (GDALGetRasterBand( poDataset, 1 ), &bGotNodata);
    if (!bGotNodata){
        logc.warn ("main", "NODATA field not available in the raster input. Using default value [0]");
        dfNoData = 0.0;
    }
    logc.debug("geotiff", "geoTIFF summary information: ****************************************");
    cout << "\tCanvas size:     \t[" << xSize << " x " << ySize << "] = [" << yellow << xSize * ySize << reset << "] pixels" << endl;
    cout << "\tPixel resolution:\t[" << yellow << xResolution << " x " << yResolution << reset << "] meter / pixel" << endl;
    cout << "\t\t> Nominal area:\t[" << cyan << xSize*ySize*xResolution*yResolution << reset << "] m2" << endl;
    cout << "\tNo data:         \t[" << (bGotNodata ? green : red) << dfNoData << reset << "]" << endl;

    logc.info ("geotiff", "Reading");
    // load data into memory
    float **apData; //pull 2D float matrix containing the image data for Band 1
    apData = inputGeotiff.GetRasterBand(1);

    double acum = 0;
    int nodataCount = 0;
    double dataCount = 0;

    for (int row = 0; row < ySize; row++){
        for (int col = 0; col < xSize; col++){
            double z = apData[row][col];
            if (z != dfNoData){
                acum += z;
                dataCount++;
            }
            else{
                nodataCount++;
            }
        }
    }

    acum /= dataCount;


    cout << "Total of no_data:\t" << nodataCount << endl;
    cout << "Total of data:\t" << dataCount << endl;
    cout << "Mean value:\t" << acum << endl;

	return 0;


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
