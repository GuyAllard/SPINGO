/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: classifier.h
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


#ifndef __CLASSIFIER_H__
#define __CLASSIFIER_H__

#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>

#include "database.h"
#include "fasta.h"
#include "common.h"

// configuration options
struct ClassifierOptions {
    int kmerSize;
    int numThreads;
    int numBootstrap;
    int subsample;
    std::string dbFilename;
    std::string inputFilename;
    bool saveIndex;
    bool dumpAmbiguous;
};

class Classifier
{
private:
    Database referenceData_;
    Kmerizer kmerizer_;
    unsigned int kmerSize_;
    unsigned int numRefSeqs_;
    unsigned int numThreads_;
    unsigned int numBootstrap_;
    unsigned int subsampleSize_;
    bool outputAmbiguous_;
    
    boost::mutex mutex_;
    boost::thread_group threads_;

    void runThread(FastaReader &reader);
    std::vector<float> getBootstrap(const KmerSequence &querySeq, RandomGen &generator, const searchHit& hit);

public:
    Classifier( const ClassifierOptions &options );

    void classify(const std::string &queryFileName);
};

#endif /* __CLASSIFIER_H__ */
