/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: common.h
Author: Guy Allard
Copyright (C) 2014  Department of Microbiology, University College Cork 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
=============================================================================*/

#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdexcept>

#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "boost/random.hpp"

//======================================================================
// FileOpenException class
// Raised when a file cannot be opened
//======================================================================
class FileOpenException : public std::ios_base::failure
{
public:
    FileOpenException(const std::string &message) : std::ios_base::failure(message) {}
};


//======================================================================
// ScopedTimer class
// Outputs to stderr the elapsed time (in seconds) between its
// construction and destruction
//======================================================================
class ScopedTimer
{
private:
    struct timeval startClock;

public:
    ScopedTimer() { gettimeofday(&startClock, NULL); }
    ~ScopedTimer() {
        struct timeval endClock;
        gettimeofday(&endClock, NULL);
        double execTime = (endClock.tv_sec - startClock.tv_sec);
        execTime += (endClock.tv_usec - startClock.tv_usec) / 1000000.0;
        std::stringstream s;
        s << std::setprecision(5) << "(" << execTime << "s)\n";
        std::cerr << s.str();
    }
};


//======================================================================
// RandomGen
// wrapper for the boost mt19937 rng which can be used by std::shuffle
//======================================================================
struct RandomGen
{
    boost::mt19937 &generator_;
    unsigned operator()(unsigned i)
    {
        boost::uniform_int<> rng(0, i-1);
        return rng(generator_);
    }

    RandomGen(boost::mt19937 &generator) : generator_(generator) {}
};


//======================================================================
// older boost versions require this for serialization
//======================================================================
template<class archive_type, class temporal_type>
void save_to(archive_type& ar, const temporal_type& tt)
{
    ar << tt;
}


#endif /* __COMMON_H__ */
