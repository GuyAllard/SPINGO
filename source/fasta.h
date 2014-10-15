/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: fasta.h
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

#ifndef __FASTA_H__
#define __FASTA_H__

#include <string>
#include <fstream>
#include <boost/thread/mutex.hpp>

/* Fasta file handling */

//======================================================================
// DnaSequence struct
// holds the header and sequence string
//======================================================================
struct DnaSequence
{
    std::string header;
    std::string sequence;

    operator bool() const { return !header.empty() && !sequence.empty(); }
};

//======================================================================
// FastaReader class
// Reads sequences one at a time from the specified fasta file
//======================================================================
class FastaReader
{
private:
    typedef boost::mutex::scoped_lock scoped_lock;
    boost::mutex mutex_;
    std::ifstream input_;
    DnaSequence curSeq_;
    long numRead_;

public:
    FastaReader(const std::string& filename);
    ~FastaReader();
    DnaSequence readSequence();
    long numRead();
};

#endif /* __FASTA_H__ */
