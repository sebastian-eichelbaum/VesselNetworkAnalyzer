//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#ifndef NOGO_RANGE_H
#define NOGO_RANGE_H

#include <iterator>

namespace nogo
{
    /**
     * Class to implement Python-like ranges in C++. Based on the code found here:
     * https://xanduchene.wordpress.com/2013/03/17/pythonic-ranges-in-c11/
     */
    class range
    {
    private:
        const int rbegin;
        const int rend;
        int step_end;
        const int step;

    public:
        range(int begin, int end, int stepSize = 1): rbegin(begin), rend(end), step(stepSize)
        {
            if ((rend - rbegin) % step == 0)
            {
                step_end = end;
            }
            else
            {
                int nsteps = static_cast< int >((rend - rbegin) / step);
                step_end = step * (nsteps + 1) + begin;
            }
        }

        explicit range(int end): range(0, end)
        {
        }

        class iterator: public std::iterator< std::random_access_iterator_tag, int >
        {
        private:
            int c;
            range& parent;

        public:
            iterator(int start, range& p): c(start), parent(p)
            {
            }

            int operator*()
            {
                return c;
            }

            const iterator* operator++()
            {
                c += parent.step;
                return this;
            }

            iterator operator++(int)
            {
                c += parent.step;
                return iterator(c - parent.step, parent);
            }

            bool operator==(const iterator& other)
            {
                return c == other.c;
            }

            bool operator!=(const iterator& other)
            {
                return c != other.c;
            }

            iterator operator+(int s)
            {
                return iterator(parent.step * s + c, parent);
            }

            iterator operator-(int s)
            {
                return iterator(c - parent.step * s, parent);
            }

            const iterator* operator--()
            {
                c -= parent.step;
                return this;
            }

            iterator operator--(int)
            {
                c -= parent.step;
                return iterator(c - parent.step, parent);
            }
        };

        iterator begin()
        {
            return iterator(rbegin, *this);
        }
        iterator end()
        {
            return iterator(step_end, *this);
        }

        int operator[](int s)
        {
            return rbegin + s * step;
        }

        int size()
        {
            return static_cast< int >((rend - rbegin) / step);
        }
    };
} // namespace nogo

#endif // NOGO_RANGE_H
