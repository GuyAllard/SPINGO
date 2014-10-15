/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: fasta.cpp
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

#include "fasta.h"
#include "common.h"

// construct and open the fasta file
FastaReader::FastaReader(const std::string& filename) : numRead_(0)
{
    curSeq_.header = curSeq_.sequence = "";
    numRead_ = 0;
    input_.open(filename.c_str());
    if(!input_.is_open())
    {
        std::ostringstream msg;
        msg << "Could not open " << filename;
        throw FileOpenException(msg.str());
    }
}


// destroy and close the fasta file
FastaReader::~FastaReader()
{
    if(input_.is_open())
        input_.close();
}

// read a sequence from the fasta file
DnaSequence FastaReader::readSequence()
{
    std::string line;
    DnaSequence seq;
    scoped_lock lock(mutex_);
    while(std::getline( input_, line ))
    {
        if (line[0] == '>')
        {
            if (!curSeq_.header.empty())
            {
                seq.header = curSeq_.header;
                seq.sequence = curSeq_.sequence;
                curSeq_.header = line.substr(1);
                curSeq_.sequence = "";
                numRead_ ++;
                return seq;
            }
            else
            {
                curSeq_.header = line.substr(1);
                curSeq_.sequence = "";
            }
        }
        else
        {
            curSeq_.sequence.append(line);
        }
    }

    if(!curSeq_.header.empty())
    {
        seq.header = curSeq_.header;
        seq.sequence = curSeq_.sequence;
        curSeq_.header = "";
        curSeq_.sequence = "";
        numRead_ ++;
        return seq;
    }

    return DnaSequence();
}


// get the number of sequences read so far
long FastaReader::numRead()
{
    scoped_lock lock(mutex_);
    return numRead_;
}


