/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: kmerizer.h
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

#ifndef __KMERIZER_H__
#define __KMERIZER_H__

#include <string>
#include <vector>
#include <stdint.h>

struct DnaSequence;

// 32-bit int should allow for a max kmer-size of 15 including Ns, but
// memory limitations may reduce this number depending on hardware
typedef uint32_t kmerSize_t; 


//======================================================================
// KmerSequence struct 
// Holds the k-mer representation of a DNA sequence
//======================================================================
struct KmerSequence
{
    std::string header;
    std::vector<kmerSize_t> kmers;
};


//======================================================================
// Kmerizer class
// Converts a DNA string into a list of unique k-mers
//======================================================================
class Kmerizer
{
private:
    kmerSize_t kmerSize_;
    kmerSize_t numKmers_;
    kmerSize_t kmerIndex( const std::string& kmer ) const;
    kmerSize_t revCompIndex( const kmerSize_t& idx ) const;

public:
    Kmerizer() : kmerSize_(0), numKmers_(0) {}
    Kmerizer(kmerSize_t size);
    
    void setKmerSize(kmerSize_t size);
    KmerSequence kmerize(const DnaSequence& sequence) const;
    kmerSize_t numKmers() const { return numKmers_; }
    kmerSize_t kmerSize() const { return kmerSize_; }

    KmerSequence revComp(const KmerSequence &kmerSeq) const;  // reverse complement the sequence
};

#endif /* __KMERIZER_H__ */
