//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#ifndef NOGO_TYPES_H
#define NOGO_TYPES_H

#include <set>
#include <vector>

// shared_ptr
#include <memory>

// NOTE: This file is mostly used for type abbreviations to make life easier.
// Thanks to the new C++11 "using" statement, you can also define "templated" type aliases. Example: alias a vector of
// shared_ptr:
//  template< typename T >
//  using SPtrVec = std::vector< std::shared_ptr< T > >
// ... Nice, isn't it?

namespace nogo
{
    /**
     * Use this type for floating point operations
     */
    using Real = float;

    /**
     * Alias for abbreviating the often used std::shared_ptr< Type >.
     *
     * \tparam T the type to embed into the shared_ptr.
     */
    template < typename T >
    using SPtr = std::shared_ptr< T >;

    /**
     * Alias for abbreviating the often used std::unique_ptr< Type >.
     *
     * \tparam T the type to embed into the shared_ptr.
     */
    template < typename T >
    using UPtr = std::unique_ptr< T >;

    /**
     * Alias for abbreviating the often used std::shared_ptr< Type >.
     *
     * \tparam T the type to embed into the shared_ptr.
     */
    template < typename T >
    using ConstSPtr = std::shared_ptr< const T >;

    /**
     * Alias for abbreviating the often used std::weak_ptr< Type >.
     *
     * \tparam T the type to embed into the weak_ptr.
     */
    template < typename T >
    using WPtr = std::weak_ptr< T >;

    /**
     * Alias for abbreviating the often used std::weak_ptr< Type >.
     *
     * \tparam T the type to embed into the weak_ptr.
     */
    template < typename T >
    using ConstWPtr = std::weak_ptr< const T >;

    /**
     * Alias for abbreviating often used shared pointer vector containers
     *
     * \tparam T the type in the shared_ptr vector.
     */
    template < typename T >
    using SPtrVec = std::vector< std::shared_ptr< T > >;

    /**
     * Alias for abbreviating often used shared pointer vector containers
     *
     * \tparam T the type in the shared_ptr vector.
     */
    template < typename T >
    using ConstSPtrVec = std::vector< std::shared_ptr< const T > >;

    /**
     * Alias for abbreviating often used shared pointer set containers
     *
     * \tparam T the type in the shared_ptr vector.
     */
    template < typename T >
    using SPtrSet = std::set< std::shared_ptr< T > >;

    /**
     * Alias for abbreviating often used shared pointer set containers
     *
     * \tparam T the type in the shared_ptr vector.
     */
    template < typename T >
    using ConstSPtrSet = std::set< std::shared_ptr< const T > >;
} // namespace nogo

#endif // NOGO_TYPES_H
