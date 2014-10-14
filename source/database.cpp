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
    numLevels_ = 0;
    index_.clear();
    
    try {
        index_.resize(kmerizer.numKmers() + 1);
    }
    catch(std::bad_alloc&) {
        std::cerr << "\nError: Could not allocate memory for kmer size " << kmerizer.kmerSize() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    annotations_.clear();
    annotationIds_.clear();
   
    FastaReader reader(fileName);
    ScopedTimer tim;
    std::cerr << "Loading reference database: " << fileName.c_str() << "\n";
    
    // read first entry separately, use it to define number of taxonomic levels
    DnaSequence seq = reader.readSequence();
    if(!seq)
    {
        std::cerr << "\nError: incorrect database format";
        exit(EXIT_FAILURE);
    }
    numLevels_ = std::count(seq.header.begin(), seq.header.end(), '\t');

    if (numLevels_ == 0)
    {
        std::cerr << "\nError: No taxonomic levels defined";
        exit(EXIT_FAILURE);
    }

    annotationIds_.resize(numLevels_);

    addSequence(kmerizer.kmerize(seq));

    // now load up the rest of the reference sequences
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


// add a sequence to the database
void Database::addSequence(const KmerSequence & seq)
{
    boost::mutex::scoped_lock lock(mutex_);
    parseHeader(seq.header);
    
    for (std::vector<kmerSize_t>::const_iterator it=seq.kmers.begin(); it != seq.kmers.end(); ++it)
        index_[*it].push_back(numSequences_);
    
    numSequences_ ++;
}


// thread that actually does the loading
void Database::loadThread(FastaReader &reader, const Kmerizer& kmerizer)
{
    DnaSequence seq;
    while(seq = reader.readSequence())
        addSequence(kmerizer.kmerize(seq));
}


// parse the sequence header and extract species and genus information
void Database::parseHeader(const std::string& header)
{
    // check that we have the required number of fields
    // and store the tab positions
    std::vector<size_t> tabPositions;
    size_t tabPos, startPos = 0;
    while((tabPos = header.find("\t", startPos)) != std::string::npos)
    {   
        tabPositions.push_back(tabPos);
        startPos = tabPos + 1;
    }

    if (tabPositions.size() != numLevels_)
    {
        std::cerr << "\nError: Invalid sequence header - incorrect number of taxonomic levels\n" << header << std::endl;
        exit(EXIT_FAILURE);
    }

    for (unsigned int i=0; i<tabPositions.size(); ++i)
    {
        size_t startPos = tabPositions[i]+1;
        std::string annotation;
        if (i == tabPositions.size() -1)
        {
            annotation = header.substr(startPos);
        }
        else
        {
            annotation = header.substr(startPos, tabPositions[i+1] - startPos);
        }

        std::vector<std::string>::iterator findit = std::find(annotations_.begin(), annotations_.end(), annotation);
        if (findit != annotations_.end())
        {
            annotationIds_[i].push_back(std::distance(annotations_.begin(), findit));
        }
        else
        {
            annotations_.push_back(annotation);
            annotationIds_[i].push_back(annotations_.size() - 1);
        }
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
    hit.annotationIds.resize(numLevels_);

    std::vector<unsigned int>::iterator searchIt=scores.begin();

    while((searchIt = std::find(searchIt, scores.end(), topScore)) != scores.end())
    {
        int id = std::distance(scores.begin(), searchIt);
        for(unsigned int i=0; i<numLevels_; i++)
            hit.annotationIds[i].push_back(annotationIds_[i][id]);
        ++searchIt;
    }

    hit.score = static_cast<float>(topScore) / static_cast<float>(query.kmers.size());
    return hit;
}


const std::string& Database::annotationFromId(const unsigned int id) const
{
    return annotations_[id];
}
