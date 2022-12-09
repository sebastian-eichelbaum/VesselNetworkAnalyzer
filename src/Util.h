//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#ifndef NOGO_UTILS_H
#define NOGO_UTILS_H

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <iostream>

#include "Types.h"

#include "Logger.h"
#define LogTag "nogo/Util"

namespace nogo
{
    /**
     * Split a given string into tokens by using a delimiter char.
     *
     * \param theString the string to split
     * \param delim the delimiter, space by default
     *
     * \return the split string elements. Is the provided string if cannot be split.
     */
    inline std::vector< std::string > split(const std::string& theString, const char& delim = ' ')
    {
        // NOTE: this can be done nicely with std::sregex_token_iterator but crashes on GCC 4.9

        // So use a string stream
        std::stringstream ss(theString);
        std::string item;
        std::vector< std::string > elems;
        while (std::getline(ss, item, delim))
        {
            elems.push_back(item);
        }
        return elems;
    }

    /**
     * Implement basic byte order swap.
     *
     * \tparam T the type of the data
     * \param value the value to swap the order for.
     * \param swapSize the size of the swap-block -> Example: byte order 1 2 3 4 5 6 7 8 with swapSize = 4 will become 4
     * 3 2 1 8 7 6 5. Very useful for types made up of other, trivial types.
     *
     * \return the swapped byte order value.
     */
    template < typename T >
    T byteSwap(const T& value, const size_t swapSize = sizeof(T))
    {
        T res; // NOTE: calls constructor. Maybe use an in-place swap.

        // Interpret value and result as a bunch of bytes:
        char const* const in = reinterpret_cast< char const* const >(&value);
        char* const out = reinterpret_cast< char* >(&res);

        // Iterate and swap
        size_t size = sizeof(T);
        for (size_t i = 0; i < size; ++i)
        {
            out[i] = in[static_cast< int >(std::floor(i / swapSize)) * swapSize + swapSize - (i % swapSize) - 1];
        }
        return res;
    }


    /**
     * Wait for the user to press any key.
     */
    inline void awaitUser(const std::string& msg = "")
    {
        do
        {
            std::cout << msg << '\n' << "Press a key to continue...";
        } while (std::cin.get() != '\n');
    }
} // namespace nogo

#endif // NOGO_UTILS_H
