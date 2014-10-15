/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: kmerizer.cpp
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

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "kmerizer.h"
#include "fasta.h"

// construct with a given kmer size
Kmerizer::Kmerizer(kmerSize_t size)
{
    setKmerSize(size);
}

// set a new kmer size
void Kmerizer::setKmerSize(kmerSize_t size)
{
    kmerSize_ = size;
    numKmers_ = 1 << (kmerSize_ << 1);
}

// convert a DNA sequence into a kmer sequence
KmerSequence Kmerizer::kmerize(const DnaSequence& sequence) const
{
    kmerSize_t seqlen = sequence.sequence.length();
    assert(seqlen >= kmerSize_);
    assert(numKmers_ > 0);
    kmerSize_t numKmersInSeq = seqlen - kmerSize_;

    KmerSequence kmerseq;
    kmerseq.header = sequence.header;

    for(kmerSize_t i=0; i <= numKmersInSeq; ++i)
    {
        kmerSize_t idx = kmerIndex(sequence.sequence.substr(i, kmerSize_));
        if(std::find(kmerseq.kmers.begin(), kmerseq.kmers.end(), idx) == kmerseq.kmers.end())
            kmerseq.kmers.push_back(idx);
    }

    return kmerseq;
}

// convert a kmer to its associated index
kmerSize_t Kmerizer::kmerIndex( const std::string& kmer ) const
{
    assert(kmer.length() == kmerSize_);

    kmerSize_t kmerId = 0;

    for(kmerSize_t i = 0; i < kmerSize_; i++)
    {
        unsigned int i2 = i << 1;
        switch (kmer[i])
        {
            case 'A':
            case 'a':
                // kmerId |= 0 << i2; // unnecessary
                break;
            case 'C':
            case 'c':
                kmerId |= 1 << i2;
                break;
            case 'G':
            case 'g':
                kmerId |= 2 << i2;
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                kmerId |= 3 << i2;
                break;
            default:
                return numKmers_; // non ACGTU base encountered
        }
    }
    return kmerId;
}

// convert a kmer sequence into its reverse complement
KmerSequence Kmerizer::revComp(const KmerSequence &kmerSeq) const
{
    KmerSequence revComp;
    revComp.header = kmerSeq.header;
    revComp.kmers.reserve(kmerSeq.kmers.size());
    for(std::vector<kmerSize_t>::const_iterator it = kmerSeq.kmers.begin(); it != kmerSeq.kmers.end(); ++it)
        revComp.kmers.push_back(revCompIndex(*it));
    
    return revComp;
}


// reverse complement a single kmer
kmerSize_t Kmerizer::revCompIndex(const kmerSize_t &idx) const
{
   // early-out if the original kmer is an 'N' kmer
   if (idx == numKmers_)
       return idx;

   kmerSize_t mask = 3;     // bitmask
   kmerSize_t comp = ~idx;  // bitwise complement
   kmerSize_t revC = 0;     // result
   
   // reverse the complemented index
   for (unsigned int i=0; i<kmerSize_; i++)
   {
       revC <<= 2;          // make space for the next bit-pair
       revC |= comp & mask; // add the bit-pair from the complement
       comp >>= 2;          // select the next bit-pair of the complement
   }

   return revC;
}
