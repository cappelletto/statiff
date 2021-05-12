/**
 * @file statiff.hpp
 * @author Jose Cappelletto (cappelletto@gmail.com)
 * @brief Base namespace definition
 * @version 0.1
 * @date 2020-10-31
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef _STATIFF_HPP_

#define _STATIFF_HPP_

namespace statiff{

    enum ErrorCode{
        NO_ERROR    =  0,    // Standard informative message
        INVALID_ARG = -1,    // non-critical warning message
        INVALID_DATA= -2     // debug (verbose) message
    };

    enum Units{
        UNIT_PIXEL  = 0,    // default, dimensions and positions are given in pixels
        UNIT_MM     = 1,    // dimensions are given in milimeters
        UNIT_CM     = 2,    // dimensions are given in centimeters
        UNIT_M      = 3,    // dimensions are given in meters
        UNIT_PCT    = 4     // dimensions are given as a proportion of the input image size
    };

};

#endif // _STATIFF_HPP_