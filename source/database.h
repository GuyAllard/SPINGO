/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: database.h
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

#ifndef __DATABASE_H__
#define __DATABASE_H__

// STL includes
#include <vector>
#include <deque>

// support files
#include "kmerizer.h"
#include "fasta.h"

// serialization support
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/deque.hpp>

// threading
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>

// increment this if the database format changes
#define DATABASE_FILE_VERSION 1

//======================================================================
// searchHit struct
// holds the result of a database search
//======================================================================
struct searchHit
{
    float score;                            // similarity score
    std::vector<std::vector<unsigned int> > annotationIds;  // annotations associated with this hit
};


//======================================================================
// Database class
// searchable database indexed by kmer
//======================================================================
class Database
{
private:
    unsigned int numSequences_;                 // total number of sequences in the database
    unsigned int numLevels_;                    // number of taxonomic levels

    void parseHeader(const std::string& header);// parse and store sequence header

    std::deque<std::vector<unsigned int> > index_; // main kmer index
    
    std::vector<std::string> annotations_;
    std::vector<std::vector<unsigned int> > annotationIds_;
   
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        if(version != DATABASE_FILE_VERSION)
            throw boost::archive::archive_exception(boost::archive::archive_exception::unsupported_version); 

        ar & numSequences_;
        ar & numLevels_;
        ar & annotations_;
        ar & annotationIds_;
        ar & index_;
    }
    
    // threading
    boost::mutex mutex_;
    boost::thread_group threads_;
    void loadThread(FastaReader& reader, const Kmerizer& kmerizer);
    void addSequence(const KmerSequence& seq);

public:
    Database() : numSequences_(0) {}
    ~Database() {}
    
    // load database from file
    bool load(const std::string& fileName, const Kmerizer& kmerizer, const unsigned int numThreads);

    // get number of sequences in the database
    const int numSequences() { return numSequences_; }

    // search the database
    searchHit search(const KmerSequence& query);

    // get annotation string from id
    const std::string& annotationFromId(const unsigned int id) const;
    const unsigned int & numLevels() const { return numLevels_; }
};

BOOST_CLASS_VERSION(Database, DATABASE_FILE_VERSION)
#endif /* __DATABASE_H__ */
