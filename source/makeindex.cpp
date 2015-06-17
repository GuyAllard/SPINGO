/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: makeindex.cpp
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

#include <boost/program_options.hpp>
#include <fstream>
#include "common.h"
#include "database.h"

// version
static const std::string versionString("Version 1.3");

// architecture
#ifdef __x86_64__
static const std::string archString("(64bit)");
#else
static const std::string archString("(32bit)");
#endif

// defaults
static const int minKmerSize = 1;
static const int maxKmerSize = 15;
static const int defaultKmerSize = 8;
static const int minThreads = 1;
static const int defaultThreads = 1;


// somewhere to store the program options
struct ProgramOptions {
    int kmerSize;
    std::string dbFilename;
    int numThreads;
};


// process the command-line args
ProgramOptions parseCommandLine(int argc, char **argv)
{
    // parse command line options
    ProgramOptions options;
    namespace po = boost::program_options;
    std::string title("SPINDEX - SPINGO index creator.\n" + versionString + " " + archString);
    po::options_description desc("Available options");
    desc.add_options()
        (
            "help,h",
            "display this message"
        )
        (
            "version,v",
            "display version information"
        )
        (
            "kmersize,k",
            po::value<int>(&options.kmerSize)->default_value(defaultKmerSize),
            "K-mer size"
        )
        (
            "database,d",
            po::value<std::string>(&options.dbFilename),
            "path to the fasta format reference sequence database"
        )
        (
            "processors,p",
            po::value<int>(&options.numThreads)->default_value(defaultThreads),
            "number of processor threads"
        );

    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, desc), vm);
        
        if(vm.count("help"))
        {
            std::cerr << title << std::endl;
            std::cerr << desc << std::endl;
            exit(EXIT_SUCCESS);
        }

        if(vm.count("version"))
        {
            std::cerr << title << std::endl;
            exit(EXIT_SUCCESS);
        }

        po::notify(vm);
        
        // check for invalid options
        if(vm.count("kmersize") && (options.kmerSize < minKmerSize || options.kmerSize > maxKmerSize))
        {
            std::ostringstream msg;
            msg << "kmersize (--kmersize, -k) = " << options.kmerSize;
            msg << ": value must be in the range [" << minKmerSize << "," << maxKmerSize << "]";
            throw po::error(msg.str());
        }
                   
        if(vm.count("numthreads") && options.numThreads < minThreads)
        {
            std::ostringstream msg;
            msg << "numthreads (--numthreads, -n) = " << options.numThreads;
            msg << ": value must be >= " << minThreads;
            throw po::error(msg.str());
        }

        if(!vm.count("database"))
            throw po::error("database not specified");
    }
    catch(po::error& e)
    {
        std::cerr << "OOPS! " << e.what() << std::endl;
        std::cerr << "Try 'spindex --help' for available options" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    return options;
}

int main(int argc, char **argv)
{
    ProgramOptions options = parseCommandLine(argc, argv);
    
    try {
        Kmerizer kizer(options.kmerSize);
        Database db;
        db.load(options.dbFilename, kizer, options.numThreads);
        
        ScopedTimer tim;
        std::stringstream s;
        s << options.dbFilename << ".idx_" << options.kmerSize;
        std::cerr << "Writing index to " << s.str() << "....";
        std::ofstream ofs(s.str().c_str());
        boost::archive::binary_oarchive oa(ofs);
        //oa << db;
        save_to(oa, db);
        std::cerr << "done. ";
    }
    catch (FileOpenException &) {
        std::cerr << "\nError: Could not read from file" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return 0;
}
