/**
 * @file options.h
 * @brief Argument parser options based on args.hxx
 * @version 1.1
 * @date 18/06/2020
 * @author Jose Cappelletto
 */

#ifndef _PROJECT_OPTIONS_H_

#define _PROJECT_OPTIONS_H_

#include "headers.h"
#include "../external/args.hxx"
#include <iostream>

using namespace statiff;

args::ArgumentParser argParser("","");
args::HelpFlag 	     argHelp(argParser, "help", "Display this help menu", {'h', "help"});
args::CompletionFlag completion(argParser, {"complete"});	//TODO: figure out why is missing in current version of args.hxx

args::ValueFlag <std::string> 	argInput(argParser, "input", "Input bathymetry map. TIFF file or XYZ point collection", {"input"});
// args::Positional<std::string> 	argInput(argParser,     "input",    "Input bathymetry map. TIFF file or XYZ point collection");
args::ValueFlag	<std::string> 	argOutput(argParser,    "output",   "Output file",{'o',"output"});
args::ValueFlag	<int> 	        argVerbose(argParser,   "verbose",  "Define verbosity level",                                                   {"verbose"});
// Free parameters for debugging
args::ValueFlag	<int> 	argNThreads(argParser,  "nthreads", "Number of multithread workers. Used as a hint for OpenMP parallel blocks",  {"int"});
args::ValueFlag	<int> 	argIntParam(argParser,  "param",    "User defined parameter INTEGER for testing purposes",  {"int"});
args::ValueFlag	<float> argFloatParam(argParser,"param",    "User defined parameter FLOAT for testing purposes",    {"float"});


int initParser(int argc, char *argv[]){
        //*********************************************************************************
    /* PARSER section */
    std::string descriptionString =
        "statiff - small analysis module to extract several stats from geoTIFF maps \
    Compatible interface with geoTIFF bathymetry datasets via GDAL. \
    Multithread support via OpenMP";

    argParser.Description(descriptionString);
    argParser.Epilog("Author: J. Cappelletto (GitHub: @cappelletto)\n");
    argParser.Prog(argv[0]);
    argParser.helpParams.width = 120;

    try
    {
        argParser.ParseCLI(argc, argv);
    }
    catch (const args::Completion &e)
    {
        cout << e.what();   // not working, fail to throw anything
        return 0;
    }

    catch (args::Help)
    { // if argument asking for help, show this message
        cout << argParser;
        return ErrorCode::INVALID_ARG;
    }
    catch (args::ParseError e)
    { //if some error ocurr while parsing, show summary
        std::cerr << e.what() << std::endl;
        std::cerr << "Use -h, --help command to see usage" << std::endl;
        return ErrorCode::INVALID_ARG;
    }
    catch (args::ValidationError e)
    { // if some error at argument validation, show
        std::cerr << "Bad input commands" << std::endl;
        std::cerr << "Use -h, --help command to see usage" << std::endl;
        return ErrorCode::INVALID_ARG;
    }
    // Start parsing mandatory arguments
    if (!argInput)
    {
        cerr << "Mandatory <input> file name missing" << endl;
        cerr << "Use -h, --help command to see usage" << endl;
        return ErrorCode::INVALID_ARG;
    }

    cout << cyan << "statiff" << reset << endl; // CREATE OUTPUT TEMPLATE STRING
    cout << "\tOpenCV version:\t" << yellow << CV_VERSION << reset << endl;
    cout << "\tGit commit:\t" << yellow << GIT_COMMIT << reset << endl
         << endl;
    cout << "\tBuilt:\t" << __DATE__ << " - " << __TIME__ << endl;   // TODO: solve, make is complaining about this
    return 0;
}

#endif //_PROJECT_OPTIONS_H_