/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: database.cpp
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

#include "database.h"
#include "fasta.h"
#include "common.h"
#include <algorithm>
#include <iostream>
#include <math.h>

// load the reference sequences and convert to kmer index
bool Database::load(const std::string& fileName, const Kmerizer& kmerizer, const unsigned int numThreads)
{
    numSequences_ = 0;
    index_.clear();
    
    try {
        index_.resize(kmerizer.numKmers() + 1);
    }
    catch(std::bad_alloc&) {
        std::cerr << "\nError: Could not allocate memory for kmer size " << kmerizer.kmerSize() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    species_.clear();
    speciesIds_.clear();
    genera_.clear();
    genusIds_.clear();
    
    FastaReader reader(fileName);
    ScopedTimer tim;
    std::cerr << "Loading reference database: " << fileName.c_str() << "\n";
    
    for (unsigned int i=0; i<numThreads; i++)
        threads_.create_thread(boost::bind(&Database::loadThread, this, boost::ref(reader), boost::ref(kmerizer)));

    threads_.join_all();
  
    // something wrong? 
    if(numSequences_ == 0)
    {
        std::ostringstream msg;
        msg << "Could not read from database " << fileName;
        throw FileOpenException(msg.str());
    }
    
    // remove the index for k-mers containing non-ATGCU bases
    index_[kmerizer.numKmers()].clear();

    // all done!
    std::cerr << numSequences_ << " reference sequences loaded ";
    return true;
}


// thread that actually does the loading
void Database::loadThread(FastaReader &reader, const Kmerizer& kmerizer)
{
    DnaSequence seq;
    while(seq = reader.readSequence())
    {
        KmerSequence kmerSeq = kmerizer.kmerize(seq);
        
        boost::mutex::scoped_lock lock(mutex_);
        parseHeader(kmerSeq.header);
        
        for (std::vector<kmerSize_t>::const_iterator it=kmerSeq.kmers.begin(); it != kmerSeq.kmers.end(); ++it)
        {
            index_[*it].push_back(numSequences_);
        }
        numSequences_ ++;
    }
}


// parse the sequence header and extract species and genus information
void Database::parseHeader(const std::string& header)
{
    size_t sppTabPos = header.find("\t");
    size_t genusTabPos = header.find("\t", sppTabPos + 1);

    if ((sppTabPos == std::string::npos) || (genusTabPos == std::string::npos))
    {
        std::cerr << "\nError: Invalid sequence header -\n" << header << std::endl;
        exit(EXIT_FAILURE);
    }
    
    unsigned int sppStart = sppTabPos + 1;
    unsigned int genusStart = genusTabPos + 1;

    std::string species = header.substr(sppStart, genusTabPos - sppStart);
    std::string genus = header.substr(genusStart);
   
    std::vector<std::string>::iterator findit = std::find(species_.begin(), species_.end(), species);
    if (findit != species_.end())
        speciesIds_.push_back(std::distance(species_.begin(), findit));
    else
    {
        species_.push_back(species);
        speciesIds_.push_back(species_.size() - 1);
    }
    
    findit = std::find(genera_.begin(), genera_.end(), genus);
    if (findit != genera_.end())
        genusIds_.push_back(std::distance(genera_.begin(), findit));
    else
    {
        genera_.push_back(genus);
        genusIds_.push_back(genera_.size() - 1);
    }
}


searchHit Database::search(const KmerSequence& query)
{
    std::vector<unsigned int> scores(numSequences_, 0);
    
    unsigned int topScore = 0;
    for (std::vector<kmerSize_t>::const_iterator it=query.kmers.begin(); it != query.kmers.end(); ++it)
    {
        std::vector<unsigned int>* seqList = &index_[*it];
        for (std::vector<unsigned int>::iterator it2 = seqList->begin(); it2 != seqList->end(); ++it2)
        {
            if(++scores[*it2] > topScore)
                topScore = scores[*it2];
        }
    }

    searchHit hit;
    std::vector<unsigned int>::iterator searchIt=scores.begin();
    while((searchIt = std::find(searchIt, scores.end(), topScore)) != scores.end())
    {
        int id = std::distance(scores.begin(), searchIt);
        hit.genusIds.push_back(genusIds_[id]);
        hit.speciesIds.push_back(speciesIds_[id]);
        ++searchIt;
    }

    hit.score = static_cast<float>(topScore) / static_cast<float>(query.kmers.size());
    return hit;
}

// get the species string from a species index value
const std::string& Database::speciesFromId(const unsigned int id) const
{
    return species_[id];
}

// get the genus string from a genus index value
const std::string& Database::genusFromId(const unsigned int id) const
{
    return genera_[id];
}

