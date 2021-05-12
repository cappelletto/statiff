/**
 * @file    statiff.cpp
 * @author  Jose Cappelletto (cappelletto@gmail.com)
 * @brief   geoTIFF statistics calculator tool. Use to analyze raster geoTIFF generated during the calculation of Landability and Measurability scores
 *          gdalinfo does not provide control over the histogram and the statistics calculation
 * @version 0.2
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
using namespace chrono;
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
    string  inputFileName    = ""; // command arg defined
    string outputFileName    = ""; // command arg defined
    bool outputConsoleOnly = false;

    Units roiUnits = UNIT_PIXEL;
    std::string unitString = "pixels";
    int verbosity = 0;
    if (argVerbose) verbosity = args::get(argVerbose);

    const static std::unordered_map<std::string,int> unit_map{
        {"px",  UNIT_PIXEL},
        {"mm",  UNIT_MM},
        {"cm",  UNIT_CM},
        {"m",   UNIT_M},
        {"pc",  UNIT_PCT},
    };

    if (argUnits){
        // parsing user defined units
        try{
            switch(unit_map.at(args::get(argUnits))){
                case UNIT_PIXEL:
                    roiUnits = UNIT_PIXEL;
                    unitString = "pixels";
                    break;
                case UNIT_MM:
                    roiUnits = UNIT_MM;
                    unitString = "mm";
                    break;
                case UNIT_CM:
                    roiUnits = UNIT_CM;
                    unitString = "cm";
                    break;
                case UNIT_M:
                    roiUnits = UNIT_M;
                    unitString = "m";
                    break;
                case UNIT_PCT:
                    roiUnits = UNIT_PCT;
                    unitString = "percent";
                    break;
                default:
                    roiUnits = UNIT_PIXEL;
                    logc.warn("unit", "User defined unit not recognized. Using default [pixels]");
                    break;
            }
        }
        catch (exception &ex){
            roiUnits = UNIT_PIXEL;
            logc.warn("unit", "User defined unit not recognized. Using default [pixels]");
        }
        if (verbosity> 0){
            s << "User defined unit: [" << yellow << unitString << reset << "]";
            logc.info("config", s);
        }
    }
    else{
        roiUnits = UNIT_PIXEL;
        if (verbosity> 0)   logc.info("config", "Using default dimension units [pixels]");
    }

    if (argInput) inputFileName = args::get(argInput); //input file is mandatory argument.
    if (inputFileName.empty()){ //not defined as command line argument? let's use config.yaml definition
        // ERROR! We do not have any definition of the input file
        logc.error ("main", "Input file missing. Please define it using --input='filename'");
        return ErrorCode::INVALID_ARG;
    }

    if (argOutput) outputFileName = args::get(argOutput); //input file is mandatory argument.
    // TODO: IF NO OUTPUT FILE IS PROVIDED, WE CAN ASSUME THAT WE ARE ASKING FOR CONSOLE OUTPUT ONLY (BATCH PROCESS)
    if (outputFileName.empty()){ //not defined as command line argument? let's use config.yaml definition
        // ERROR! We do not have any definition of the input file
        if (verbosity >= 1){
            logc.info ("main", "Output file missing, console output only. If you need to export the results please define an output file using [--output]");
        }
        outputConsoleOnly = true; // no outputfile provided, exporting results to console
        // return ErrorCode::INVALID_ARG;
    }

    double histMin = 0.0, histMax = 1.0;
    int nBins = 100;
    if (!argNoHistogram)
    {
        if (argHistMax)   histMax = args::get(argHistMax);
        if (argHistMin)   histMin = args::get(argHistMin);
        if (argHistBin)     nBins = args::get(argHistBin);
        // if argNoHistogram is called, these paremeters will be ignored

        // validating input:
        if (histMin >= histMax){
            s << "Wrong histogram range arguments. histMin [" << yellow << histMin << reset << "] must be lower than histMax [" << yellow << histMax << reset << "]";
            logc.error("main", s);
            return statiff::INVALID_ARG;
        }
        
        if (nBins <= 0){
            s << "Number of bins nBins[" << yellow << nBins << reset << "] must be positive integer";
            logc.error("main", s);
            return statiff::INVALID_ARG;
        }
    }

    // CHECK IF ROI DEFINED
    bool bRoiDefined = false;
    long int roiWidth = -1, roiHeight = -1;
    if (argRoiHeight){
        bRoiDefined = true;
        roiHeight = args::get(argRoiHeight);
    }
    if (argRoiWidth){
        bRoiDefined = true;
        roiWidth = args::get(argRoiWidth);
    }

    // PRINT SUMMARY
    if (verbosity >= 1){
        logc.debug("main", "Summary information *****************************************************");
        cout << "\tInput file:   \t" << green << inputFileName << reset << endl;
        cout << "\tOutput file:  \t" << green << outputFileName << reset << endl;
        cout << "\tExport header:\t" << yellow << (argNoHeader ? "true" :  "false") << reset << endl;
        if (argNoHistogram)
            cout << yellow << "No histogram will be computed (--nohist option passed)" << reset << endl;
        else
            cout << "\tHistogram:    \t" << green << "[" << histMin << " : " << histMax << "]" << reset << "\t #Bins: " << yellow << nBins << reset << endl;
        if (bRoiDefined){
            cout << "ROI defined: [" << yellow << roiWidth << " x " << roiHeight << reset << "] " << unitString << endl;
        }
    }

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
    int *dimensions;
    dimensions = inputGeotiff.GetDimensions();

    // canvas dimension
    long int xSize = dimensions[0]; //NCols
    long int ySize = dimensions[1]; //nRows
    // long int ySize = poDataset->GetRasterYSize();
    // pixel resolution
    double xResolution, yResolution;
	double adfGeoTransform[6]; // 6-DOF geoTIFF params: center, resolution, rotation...
	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
	{
	    xResolution = fabs(adfGeoTransform[1]);
	    yResolution = fabs(adfGeoTransform[5]);
	}
    // raster NO_DATA value
    int bGotNodata = 0;
    double dfNoData = GDALGetRasterNoDataValue (GDALGetRasterBand( poDataset, 1 ), &bGotNodata);
    if (!bGotNodata && (verbosity >=2 )){
        logc.warn ("main", "NODATA field not available in the raster input. Using as default value [-9999]");
        dfNoData = -9999.0;
    }
    double nom_area = xSize*ySize*xResolution*yResolution;
    if (verbosity >= 1){
        logc.debug("geotiff", "geoTIFF summary information: *****************************************");
        cout << "\tCanvas size:     \t[" << xSize << " x " << ySize << "] = [" << yellow << xSize * ySize << reset << "] pixels" << endl;
        cout << "\tPixel resolution:\t[" << yellow << xResolution << " x " << yResolution << reset << "] meter / pixel" << endl;
        cout << "\t\t> Nominal area:\t[" << cyan << nom_area << reset << "] m2\t[" << yellow << nom_area/10000 << reset << "] ha" << endl;
        cout << "\tNo data:         \t[" << (bGotNodata ? green : red) << dfNoData << reset << "]" << endl;
        logc.debug ("geotiff", "Loading geoTIFF data into memory. This may take a while ...");
    }
    //let's benchmark reading time
    auto chrono_start = high_resolution_clock::now(); 
    // load data into memory
    float **apData; //pull 2D float matrix containing the image data for Band 1
    apData = inputGeotiff.GetRasterBand(1);
    auto chrono_stop     = high_resolution_clock::now(); 
    auto chrono_duration = duration_cast<milliseconds>(chrono_stop - chrono_start); 
    if (verbosity >= 1){
        s << "Ellapsed time reading geoTIFF: " << blue << (chrono_duration.count()/1000.0) << reset << " seconds...";
        logc.info ("time", s);
    }    
    // Next step: compute the min, max, mean
    double acum = 0.0, acum_nz = 0.0; // _nz vars correspond to the same variables without the ZERO elements (assuming unimodal distribution)
    int nodataCount = 0;
    int cont_nz = 0;
    double _min, _max; // we cannot retrieve them from the data matrix because we need to remove first the nodata fields
    double _min_nz, _max_nz; // we cannot retrieve them from the data matrix because we need to remove first the nodata fields
    double delta = (histMax - histMin)/(double) nBins; // width of each histogram bin
    vector<double> v; //let's preallocate half of the required memory to avoid triggered realloc during exec
    v.reserve(xSize*ySize/2);
    vector<int> histogram(nBins, 0);    // histogram container

    // double z;    // too much overhead multithreading, not worth creating too many workers
    // #pragma omp for reduction(+:acum)
    chrono_start = high_resolution_clock::now();
    for (int row = 0; row < ySize; row++){
        for (int col = 0; col < xSize; col++){
            double z = apData[row][col];
            // z = apData[row][col];
            if (z != dfNoData && !isnan(z)){
                v.push_back(z);  // push into the vector
                acum += z;  // increase accumulator (to compute mean value)
                if (z != 0) cont_nz++;
                // if z = 0, it won't affect the value of the acum
                // let's compute the corresponding bin
                if (!argNoHistogram){
                    int bin;
                    if      (z <= histMin) bin = 0;       // cap min
                    else if (z >= histMax) bin = nBins-1; // cap max
                    else{
                        bin = floor((z - histMin)/delta);
                    }
                    // cout << bin << " ";
                    histogram[bin] = histogram[bin] + 1;
                }
            }
            else{
                nodataCount++;// increase counter of nodata pixels. It must match TotalPixel - N
            }
        }
    }
    chrono_stop = high_resolution_clock::now();
    chrono_duration = duration_cast< milliseconds>(chrono_stop - chrono_start); 
    if (verbosity >= 1){
        s << "Ellapsed time pre-processing data: " << blue << chrono_duration.count() << reset << " milliseconds...";
        logc.info ("time", s);
        logc.info ("time", "Now computing stats ...");
    }


    chrono_start = high_resolution_clock::now();
    int N =  v.size();
    double mean    = acum/N;
    double mean_nz = acum/cont_nz;
    double variance = 0, variance_nz = 0;
    double k_acum = 0,   k_acum_nz = 0;
    acum = 0;
    // Now, having the mean, we can compute the stdev (first the variance)
    // Main loop: ****************************************************************
    _min = _max = v[0]; // we use the first value as seed for min/max extraction
    for (auto z:v){
        register double dist = z - mean; // distance to the mean
        register double d2 = dist*dist;
        acum += d2;      // acum to the variance
        k_acum += d2*d2; // acum for the distribution kurtosis 

        register double dist_nz = z - mean_nz; // distance to the mean
        register double d2_nz = dist_nz*dist_nz;
        acum_nz += d2;      // acum to the variance
        k_acum_nz += d2_nz*d2_nz; // acum for the distribution kurtosis 

        if (z < _min)       _min = z; //let's update the min
        else if (z > _max)  _max = z; //let's update the max
    }
    variance = acum / N; // no Bessel correction because we are using the entire population. Otherwise we should use v.size()-1
    double stdev = sqrt(variance); 
    // kurtosis = (1/N)sum{(x-mean)^4} / (1/N)*sum{(x-mean)^2}^2..................... (x-mean)^2 = variance
    // kurtosis = (1/N) * (k_acum) / variance^2
    double kurtosis = k_acum / (variance * variance * N);
    // partial sort is o(N) linear complexity, guarantees the position of the midpoint element of the array
    std::nth_element(v.begin(), v.begin() + N/2, v.end());
    double median = v[N/2]; // the element in the middle of the partially sorted vector
    // Non-parametric skew S = (mean - median)/stdev
    // skewness is defined using the Peasron 2nd skew coefficient: p2c = 3*S
    double skew = 3*(mean - median)/stdev;
    double area = N*xResolution*yResolution;

    variance_nz = acum_nz / cont_nz; // no Bessel correction because we are using the entire population. Otherwise we should use N-1
    double stdev_nz = sqrt(variance_nz); 
    // kurtosis = (1/N)sum{(x-mean)^4} / (1/N)*sum{(x-mean)^2}^2..................... (x-mean)^2 = variance
    // kurtosis = (1/N) * (k_acum) / variance^2
    double kurtosis_nz = k_acum_nz / (variance_nz * variance_nz * cont_nz);

    // remove ZERO from the vector before sorting it
    // auto v_nz = v;
    v.erase(std::remove(v.begin(), v.end(), 0), v.end());
    // partial sort is o(N) linear complexity, guarantees the position of the midpoint element of the array
    std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
    double median_nz = v[v.size()/2]; // the element in the middle of the partially sorted vector
    // Non-parametric skew S = (mean - median)/stdev
    // skewness is defined using the Peasron 2nd skew coefficient: p2c = 3*S
    double skew_nz = 3*(mean_nz - median_nz)/stdev_nz;
    double area_nz = v.size()*xResolution*yResolution;
    double proportion = (double) N/(xSize*ySize);

    chrono_stop = high_resolution_clock::now();
    chrono_duration = duration_cast< milliseconds>(chrono_stop - chrono_start); 
    if (verbosity >= 1){
        s << "Ellapsed time computing stats: " << blue << chrono_duration.count() << reset << " milliseconds...";
        logc.info ("time", s);
    }

    // PRINT SUMMARY ON SCREEN
    if (verbosity >= 1){
        cout << "Total of no_data:\t" << nodataCount << endl;
        // cout << "Total of data:\t" << dataCount << endl;
        cout << "Vector size:\t" << N << endl;
        cout << "Proportion: \t" << proportion << endl;
        cout << "Actual area:\t" << area << endl;
        cout << "Mean value:\t" << mean << "\t\twithout ZERO" << "\t" << mean_nz << endl;
        cout << "Variance value:\t" << variance << "\twithout ZERO" << "\t" << variance_nz << endl;
        cout << "Stdev:    \t" << stdev << "\twithout ZERO" << "\t" << stdev_nz << endl;
        cout << "Kurtosis: \t" << kurtosis  << "\t\twithout ZERO" << "\t" << kurtosis_nz << endl;
        cout << "Skew:     \t" << skew  << "\t\twithout ZERO" << "\t" << skew_nz <<  endl;
        cout << reset << "Median value:\t" << median  << "\t\twithout ZERO" << "\t" << median_nz <<  endl;
        cout << "Min value:\t" << _min << endl;
        cout << "Max value:\t" << _max << endl;

        if (!argNoHistogram){
            cout << "Histogram: [";
            for (auto x:histogram)
                cout << x << " ";
            cout << "]" << endl;
        }
    }

    if (outputConsoleOnly){
        if (verbosity >= 1){
            logc.warn("main", "No summary files were exported. Use [--output] if required");
        }
        return NO_ERROR;
    }

    // DUMP TO OUTPUT FILE
    // ****************************************************************************************
    std::ofstream outFile("STAT_" + outputFileName);
    if (!argNoHeader){ // generate header
        outFile << "Filename\tNODATA\tXresolution\tYresolution\tN_Pixels\tData_Proportion\tArea_m2\tMin\tMax\tMean\tMedian\tStdev\tKurtosis\tSkew\tHist_min\tHist_max\tBins";
        // for (int i=0; i < nBins; i++)
        //     outFile << "\tbin_" << i;
        outFile << endl;
    }
    outFile << outputFileName << "\t" << dfNoData << "\t" << xResolution << "\t" << yResolution << "\t" << N << "\t" <<  proportion<< "\t" << area;
    outFile << "\t" << _min << "\t" << _max << "\t" << mean << "\t" << median << "\t" << stdev << "\t" << kurtosis << "\t" << skew << "\t" << histMin << "\t" << histMax << "\t" << nBins;
    // for (int i=0; i < nBins; i++)
    //     outFile << "\t" << histogram[i];
    outFile << endl;

    // Now we repeat the row-entry for Non-Zero stats
    outFile << "ZEROES_REMOVED" << "\t" << dfNoData << "\t" << xResolution << "\t" << yResolution << "\t" << v.size() << "\t" << proportion << "\t" << area_nz;
    outFile << "\t" << _min << "\t" << _max << "\t" << mean_nz << "\t" << median_nz << "\t" << stdev_nz << "\t" << kurtosis_nz << "\t" << skew_nz << "\t" << histMin << "\t" << histMax << "\t" << nBins;
    // for (int i=0; i < nBins; i++){
    //     if (i==0) histogram[i] = 0;     // we erase all ZERO entries in the histogram before exporting it
    //     outFile << "\t" << histogram[i];
    // }
    outFile << endl;

    // ****************************************************************************************
    // -- The actual histograms is dumped onto a separate file (filename_suffix, suffix: _HIST)
    if (!argNoHistogram){
    std::ofstream outFileHist("HIST_" + outputFileName);
        std::ofstream outFileBins("BINS_" + outputFileName); // file that contains the raw bin count in a column-wise fashion. Easier to merge using bash
        // first c
        if (!argNoHeader){ // generate header for histogram data with some context (linked stat file, number of bins, histogram range)
            outFileHist << "#Stat_file:\t" << outputFileName << endl;
            outFileHist << "#nBins:\t"     << nBins   << endl;
            outFileHist << "#histMin:\t"   << histMin << endl;
            outFileHist << "#histMax:\t"   << histMax << endl; 
        }
        // add X column with the center of each bin
        outFileHist << "xbin\thist\thist_nozero\thist_norm\thist_nozero_norm" << endl;
        for (int i=0; i < nBins; i++){
            register double xbin = histMin + delta/2 + delta*i; // centered by delta/2 
            outFileHist << xbin << "\t";            // central position of each bin. calculate to start from histMin + delta/2
            outFileHist << histogram[i] << "\t";            // original histogram
            outFileHist << (i ? histogram[i] : 0) << "\t";            // original histogram
            outFileHist << double(histogram[i]/(double)N) << "\t";            // original histogram
            outFileHist << (i ? double(histogram[i]/(double)v.size()) : 0.0) << endl;            // original histogram

            outFileBins << histogram[i];
            if (i<(nBins-1)) outFileBins << "\t";
                else         outFileBins << endl;
        }
    }
    //------------>END
    return ErrorCode::NO_ERROR; // everything went smooth
}
