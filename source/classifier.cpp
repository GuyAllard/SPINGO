/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: classifier.cpp
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
#include "classifier.h"

// construct the classifier
Classifier::Classifier(const ClassifierOptions &options )
{
    kmerSize_ = options.kmerSize;
    numThreads_ = options.numThreads;
    numBootstrap_ = options.numBootstrap;
    subsampleSize_ = options.subsample;
    kmerizer_.setKmerSize(kmerSize_);
    outputAmbiguous_ = options.dumpAmbiguous;
    
    std::stringstream s;
    s << options.dbFilename << ".idx_" << kmerSize_;
    try {
        // load cached index
        std::ifstream ifs(s.str().c_str());
        boost::archive::binary_iarchive ia(ifs);
        ScopedTimer tim;
        std::cerr << "Loading cached index from " << s.str();
        ia >> referenceData_;
        std::cerr << " done. ";
    }
    catch (boost::archive::archive_exception &) {
        // create index from raw sequences
        referenceData_.load(options.dbFilename, kmerizer_, numThreads_);
        
        // save index for future use
        if (options.saveIndex)
        {
            std::ofstream ofs(s.str().c_str());
            boost::archive::binary_oarchive oa(ofs);
            ScopedTimer tim;
            std::cerr << "Writing cached index to " << s.str();
            save_to(oa, referenceData_);
            std::cerr << " done. ";
        }
    }
    
    numRefSeqs_ = referenceData_.numSequences();
}

// perform classification in separate threads
void Classifier::runThread(FastaReader &reader)
{
    DnaSequence seq;
    boost::mt19937 randGen(0);
    RandomGen generator(randGen);

    while(seq = reader.readSequence())
    {
        // reproducible randomness
        randGen.seed(reader.numRead());

        // get a query sequence and convert to kmers
        KmerSequence fwdSeq = kmerizer_.kmerize(seq);
        
        // and the reverse complement
        KmerSequence revSeq = kmerizer_.revComp(fwdSeq);
        
        // search against database using both forward and reverse sequences
        searchHit fwdHit = referenceData_.search(fwdSeq);
        searchHit revHit = referenceData_.search(revSeq);

        // use the direction which gave the highest scoring hit
        KmerSequence &querySeq = fwdHit.score > revHit.score ? fwdSeq : revSeq;
        searchHit &hit = fwdHit.score > revHit.score ? fwdHit : revHit;

        std::vector<float> bootstraps(referenceData_.numLevels(), 0.f);
        std::vector<std::string> annotations(referenceData_.numLevels(), "AMBIGUOUS");
       
        // make each level annotations into a list of uniques
        for (unsigned int i=0; i< hit.annotationIds.size(); i++)
        {
            std::sort(hit.annotationIds[i].begin(), hit.annotationIds[i].end());
            std::vector<unsigned int>::iterator it;
            it = std::unique(hit.annotationIds[i].begin(), hit.annotationIds[i].end());
            hit.annotationIds[i].resize(std::distance(hit.annotationIds[i].begin(), it));
            if (hit.annotationIds[i].size() == 1)
                annotations[i] = referenceData_.annotationFromId(*hit.annotationIds[i].begin());
        }

        // bootstrap
        if(numBootstrap_ > 0)
            bootstraps = getBootstrap(querySeq, generator, hit);
       
        // output        
        std::stringstream s;
        s << std::setprecision(2) << std::fixed;
        s << querySeq.header.substr(0, querySeq.header.find("\t")) << "\t";
        s << hit.score << "\t";

        //for(unsigned int i=0; i<bootstraps.size(); i++)
        unsigned int i=bootstraps.size();
        while(i--)
        {
            s << annotations[i] << "\t";
            s << bootstraps[i];
            if(i>0)
            {
                s <<  "\t";
            }
            else
            {
                // dump ambiguous species
                if (hit.annotationIds[i].size() > 1 && outputAmbiguous_)
                { 
                    s << "\t";
                    for(std::vector<unsigned int>::iterator it = hit.annotationIds[i].begin(); it != hit.annotationIds[i].end(); ++it)
                    {
                        if (it != hit.annotationIds[i].begin())
                            s << ",";
            
                        s << referenceData_.annotationFromId(*it);
                    }
                }
            
                s << std::endl;
            }
        }
      
        boost::mutex::scoped_lock lock(mutex_);
        std::cout << s.str();
    }
}

// bootstrapping
std::vector<float> Classifier::getBootstrap(const KmerSequence& querySeq, RandomGen &generator, const searchHit& hit)
{
    KmerSequence bootstrap;
    std::vector<kmerSize_t>kmers = querySeq.kmers;
    std::vector<float> counts(referenceData_.numLevels(), 0.f);
    int bootstrap_size = kmers.size() / subsampleSize_;
    
    for (unsigned int bs=0; bs<numBootstrap_; bs++)
    {
        std::random_shuffle(kmers.begin(), kmers.end(), generator);
        bootstrap.kmers = std::vector<kmerSize_t>(kmers.begin(), kmers.begin() + bootstrap_size);
        searchHit bsHit = referenceData_.search(bootstrap);
        
        for (unsigned int i=0; i<referenceData_.numLevels(); i++)
        {
            if (hit.annotationIds[i].size() == 1)
            {
                if(std::find(bsHit.annotationIds[i].begin(), bsHit.annotationIds[i].end(), *hit.annotationIds[i].begin()) != bsHit.annotationIds[i].end())
                {
                    counts[i] += static_cast<float>(std::count(bsHit.annotationIds[i].begin(), bsHit.annotationIds[i].end(), *hit.annotationIds[i].begin())) / static_cast<float>(bsHit.annotationIds[i].size());
                }
            }
        }
    }

    std::vector<float> retVec;
    for (unsigned int i=0; i<referenceData_.numLevels(); i++)
    {
        retVec.push_back(counts[i] / static_cast<float>(numBootstrap_));
    }

    return retVec;
}


// classify sequences from queryFileName
void Classifier::classify(const std::string& queryFileName)
{
    ScopedTimer tim;
    std::cerr << "Classifying sequences....\n";
  
    FastaReader reader(queryFileName);
    for (unsigned int i=0; i<numThreads_; i++)
        threads_.create_thread(boost::bind(&Classifier::runThread, this, boost::ref(reader)));

    threads_.join_all();

    std::cerr << reader.numRead() << " sequences processed.";
}
