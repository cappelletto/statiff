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

};

#endif // _STATIFF_HPP_