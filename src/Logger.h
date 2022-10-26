//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

// A bit crude, but ensures the tag is undefined every time someone includes the logger.
#ifdef LogTag
    #undef LogTag
#endif

#ifndef NOGO_LOGGER_H
    #define NOGO_LOGGER_H

    #include <iomanip>
    #include <iostream>

// This file contains the preprocessor based logger. This is a rather simple logger, but as it uses a simple ostream
// interface, it can be replaced later without the need for changes in the code.

    #ifndef LogStream
        #define LogStream                                                                                              \
            std::cout << std::resetiosflags(std::ios_base::basefield | std::ios_base::floatfield |                     \
                                            std::ios_base::adjustfield)
    #endif

    #ifndef LogEnd
        #define LogEnd std::endl
    #endif

    #ifndef LogContinue
        #define LogContinue LogStream
    #endif

    #ifndef LogD
        #define LogD LogStream << "DEBUG [" << LogTag << "]: "
    #endif

    #ifndef LogI
        #define LogI LogStream << "INFO  [" << LogTag << "]: "
    #endif

    #ifndef LogW
        #define LogW LogStream << "WARN  [" << LogTag << "]: "
    #endif

    #ifndef LogE
        #define LogE LogStream << "ERROR [" << LogTag << "]: "
    #endif

    #ifndef LogGL
        #define LogGL LogStream
    #endif

#endif // NOGO_LOGGER_H
