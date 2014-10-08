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
Classifier::Classifier(const std::string& dbFileName,
                       const unsigned int kmerSize,
                       const unsigned int numThreads,
                       const unsigned int numBootstrap,
                       const unsigned int subsampleSize,
                       const bool saveIndex)
{
    kmerSize_ = kmerSize;
    numThreads_ = numThreads;
    numBootstrap_ = numBootstrap;
    subsampleSize_ = subsampleSize;
    kmerizer_.setKmerSize(kmerSize);
    
    std::stringstream s;
    s << dbFileName << ".idx_" << kmerSize;
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
        referenceData_.load(dbFileName, kmerizer_, numThreads_);
        
        // save index for future use
        if (saveIndex)
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
        // get a query sequence and convert to kmers
        randGen.seed(reader.numRead()); // ensure reproducible randomness
        KmerSequence fwdSeq = kmerizer_.kmerize(seq);
        
        // and the reverse complement
        seq.revComp();
        KmerSequence revSeq = kmerizer_.kmerize(seq);
        
        // search against database using both forward and reverse sequences
        searchHit fwdHit = referenceData_.search(fwdSeq);
        searchHit revHit = referenceData_.search(revSeq);

        // use the direction which gave the highest scoring hit
        KmerSequence &querySeq = fwdHit.score > revHit.score ? fwdSeq : revSeq;
        searchHit &hit = fwdHit.score > revHit.score ? fwdHit : revHit;

        std::pair<float, float> bootstrap = std::make_pair(0.f, 0.f);
        std::string species("AMBIGUOUS");
        std::string genus("AMBIGUOUS");
        
        // make unique list of hit genera
        std::sort(hit.genusIds.begin(), hit.genusIds.end());
        std::vector<unsigned int>::iterator it;
        it = std::unique(hit.genusIds.begin(), hit.genusIds.end());
        hit.genusIds.resize(std::distance(hit.genusIds.begin(), it));

        if (hit.genusIds.size() == 1)
        {
            // unique list of hit species
            std::sort(hit.speciesIds.begin(), hit.speciesIds.end());
            it = std::unique(hit.speciesIds.begin(), hit.speciesIds.end());
            hit.speciesIds.resize(std::distance(hit.speciesIds.begin(), it));
            
            // bootstrapping
            if (numBootstrap_ > 0)
                bootstrap = getBootstrap(querySeq, generator, hit);            
 
            // get genus and species names           
            genus = referenceData_.genusFromId(*hit.genusIds.begin());
            if(hit.speciesIds.size() == 1)
                species = referenceData_.speciesFromId(*hit.speciesIds.begin());
            else
                bootstrap.second=0.f;
       }
        
        // output        
        std::stringstream s;
        s << std::setprecision(2) << std::fixed;
        s << querySeq.header.substr(0, querySeq.header.find("\t")) << "\t";
        s << hit.score << "\t";
        s << genus << "\t";
        s << bootstrap.first << "\t";
        s << species << "\t";
        s << bootstrap.second << std::endl;
       
        boost::mutex::scoped_lock lock(mutex_);
        std::cout << s.str();
    }
}

// bootstrapping
std::pair<float, float> Classifier::getBootstrap(const KmerSequence& querySeq, RandomGen &generator, const searchHit& hit)
{
    KmerSequence bootstrap;
    std::vector<kmerSize_t>kmers = querySeq.kmers;
    
    float genusCount = 0;
    float speciesCount = 0;
    
    int bootstrap_size = kmers.size() / subsampleSize_;
    
    for (unsigned int i=0; i<numBootstrap_; i++)
    {
        std::random_shuffle(kmers.begin(), kmers.end(), generator);
        bootstrap.kmers = std::vector<kmerSize_t>(kmers.begin(), kmers.begin() + bootstrap_size);
        
        searchHit bsHit = referenceData_.search(bootstrap);
        
        if (std::find(bsHit.genusIds.begin(), bsHit.genusIds.end(), *hit.genusIds.begin()) != bsHit.genusIds.end())
        {
            genusCount += static_cast<float>(std::count(bsHit.genusIds.begin(), bsHit.genusIds.end(), *hit.genusIds.begin())) / static_cast<float>(bsHit.genusIds.size());
            if(std::find(bsHit.speciesIds.begin(), bsHit.speciesIds.end(), *hit.speciesIds.begin()) != bsHit.speciesIds.end())
                speciesCount += static_cast<float>(std::count(bsHit.speciesIds.begin(), bsHit.speciesIds.end(), *hit.speciesIds.begin())) / static_cast<float>(bsHit.speciesIds.size());
        }
   }

   return std::make_pair(genusCount   / static_cast<float>(numBootstrap_),
                          speciesCount / static_cast<float>(numBootstrap_));
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
