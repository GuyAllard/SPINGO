/*============================================================================
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
File: main.cpp
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
#include "common.h"
#include "classifier.h"

// version
static const std::string versionString("Version 1.3");

// platform
#ifdef __x86_64__
static const std::string archString("(64bit)");
#else
static const std::string archString("(32bit)");
#endif

// defaults
static const int minKmerSize = 1;
static const int maxKmerSize = 15;
static const int defaultKmerSize = 8;
static const int minBootstrap = 0;
static const int defaultBootstrap = 10;
static const int minThreads = 1;
static const int defaultThreads = 1;
static const int minSubsample = 1;
static const bool defaultWriteIndex = false;
static const bool defaultAmbiguousOutput = false;


// process the command-line args
ClassifierOptions parseCommandLine(int argc, char **argv)
{
    ClassifierOptions options;
    namespace po = boost::program_options;
    std::string title("SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.\n" + versionString + " " + archString);
    po::options_description desc("Available options");
    desc.add_options()
        (
            "help,h",
            "display this message"
        )
        (
            "version,v",
            "get version information"
        )
        (
            "kmersize,k",
            po::value<int>(&options.kmerSize)->default_value(defaultKmerSize),
            "K-mer size"
        )
        (
            "bootstrap,b",
            po::value<int>(&options.numBootstrap)->default_value(defaultBootstrap),
            "number of bootstrap samples"
        )
        (
            "subsample,s",
            po::value<int>(&options.subsample),
            "fraction of kmers to be subsampled for bootstrapping. Default is kmersize"
        )
        (
            "processors,p",
            po::value<int>(&options.numThreads)->default_value(defaultThreads),
            "number of processor threads"
        )
        (
            "database,d",
            po::value<std::string>(&options.dbFilename),
            "path to the fasta format reference database"
        )
        (
            "input,i",
            po::value<std::string>(&options.inputFilename),
            "path to the fasta format input file"
        )
        (
            "write-index,w",
            po::value<bool>(&options.saveIndex)->zero_tokens()->default_value(defaultWriteIndex),
            "if specified, index will be written to disk"
        )
        (
            "ambiguous,a",
            po::value<bool>(&options.dumpAmbiguous)->zero_tokens()->default_value(defaultAmbiguousOutput),
            "if specified, species which lead to an ambiguous hit will be listed"
        );

    po::positional_options_description p;
    p.add("input", -1);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

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
            
        if(vm.count("bootstrap") && options.numBootstrap < minBootstrap)
        {
            std::ostringstream msg;
            msg << "bootstrap (--boostrap, -b) = " << options.numBootstrap;
            msg << ": value must be >= " << minBootstrap;
            throw po::error(msg.str());
        }
        
        if(vm.count("numthreads") && options.numThreads < minThreads)
        {
            std::ostringstream msg;
            msg << "numthreads (--numthreads, -n) = " << options.numThreads;
            msg << ": value must be >= " << minThreads;
            throw po::error(msg.str());
        }
        
        if(!vm.count("subsample"))
            options.subsample = options.kmerSize;
        else if (options.subsample < 1)
        {
            std::ostringstream msg;
            msg << "subsample (--subsample, -s) = " << options.subsample;
            msg << ": value must be >= " << minSubsample;
            throw po::error(msg.str());
        }

        if(!vm.count("database"))
            throw po::error("database not specified");

        if(!vm.count("input"))
            throw po::error("input file not specified");
    }
    catch(po::error& e)
    {
        std::cerr << "OOPS! " << e.what() << std::endl;
        std::cerr << "Try 'spingo --help' for available options" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    return options;
}

int main(int argc, char **argv)
{
    // parse command line
    ClassifierOptions options = parseCommandLine(argc, argv);
    
    std::cerr << "Using " << options.numThreads << " thread" << (options.numThreads > 1 ? "s" : "");
    std::cerr << " with kmer size " << options.kmerSize;
    std::cerr << ", " << options.numBootstrap << " bootstrap samples";
    std::cerr << " and subsample size 1/" << options.subsample << std::endl;

    // run the classifier
    try {
        Classifier classifier(options);
        classifier.classify(options.inputFilename);
    }
    catch (FileOpenException &e) {
        std::cerr << "OOPS! " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return 0;
}
