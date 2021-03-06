## PROGRAM NAME
SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons  
         (alternatively: an olde English word meaning ‘strong beer’)  
  
  
## VERSION  
This document applies to SPINGO version 1.3  
  
  
## COPYRIGHT  
Copyright (C) 2014 Department of Microbiology, University College Cork  
  
  
## CONTACT  
[Marcus Claesson](m.claesson@ucc.ie)  
  
  
## CITING SPINGO  
Allard G, Ryan FJ, Jeffery IB, Claesson MJ. SPINGO: a rapid species-classifier for microbial amplicon sequences. BMC Bioinformatics. 2015 Oct 8;16(1):324. doi: 10.1186/s12859-015-0747-1. PubMed PMID: 26450747; PubMed Central PMCID: PMC4599320.  
  
  
## LICENSE  
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
  
  
## SYSTEM REQUIREMENTS
To run SPINGO and its supporting python scripts requires:
- A modern Linux operating system.
- Python2.7 or greater.

If the included SPINGO binaries fail to run on your system, then you will need to compile them from source. This requires additional dependencies. The full procedure and requirements can be found in the section [COMPILING FROM SOURCE](#COMPILING-FROM-SOURCE) at the end of this document.
  
  
## INSTALLATION:
The SPINGO project is [hosted on GitHub](https://GitHub.com/GuyAllard/SPINGO/)  

The SPINGO source code and executables can be obtained via two different methods:
  
Either clone the repository using git:  
```
git clone https://github.com/GuyAllard/SPINGO
```
Or download and extract the zip file using the 'Download ZIP' link of the GitHub project page.  
  
No further installation is necessary. The program can be run from within the cloned repository or from the location that the zip was extracted to.  
  
This document will assume that SPINGO is located in your home directory in a folder named `/home/your_username/SPINGO`. This location will be referred to from hereon in as 'SPINGODIR'.  
If you have placed SPINGO in a different location, use that path in place of 'SPINGODIR'.  
  
  
## GENERATION OF SPECIES DATABASE
Before SPINGO can be used to classify amplicons, a species specific database needs to be prepared.  
  
At the time of release, this database is derived from the 16S data files named 'release11_2_Archaea_unaligned.fa.gz' and 'release11_2_Bacteria_unaligned.fa.gz' which are available from the [Ribosomal Database Project](rdp.cme.msu.edu).  
  
If you already have these files, place a copy of them into the `SPINGODIR/database/` directory.  
If you do not have these files, they can be downloaded automatically as part of the following step.  
  
Change to the `SPINGODIR/database` directory.  
Execute the command:  
`make`  
  
If the 16S data files are not present, you will be asked if you wish to download them. Enter 'y' at the on-screen prompts to download these files. If you choose not to download the files automatically, you will need to download them manually, place them into 'SPINGODIR/database/' and re-run 'make' to complete the process.  
  
At the end of this process, the database folder should now contain a file named `RDP_11.2.species.fa` which SPINGO will use as its reference database.  
  
  
## RUNNING SPINGO
To check that SPINGO will run on your system, change into 'SPINGODIR' and execute the command:  
`./spingo -h`  
  
If this command does not output usage information to the screen, you will need to compile the program from source. See the section [COMPILING FROM SOURCE](#COMPILING-FROM-SOURCE) at the end of this document for further instructions.  
  
To run spingo with sensible default values, on a hypothetical amplicon sequence file named reads.fa, using the species specific database created previousl, and saving the output to a file named results.out:  
`SPINGODIR/spingo -d SPINGODIR/database/RDP_11.2.species.fa -i path/to/reads.fa > path/to/results.out`  
  
Note that spingo writes to stdout (the screen), so the ‘>’ operator is used to redirect this output to a file.  
  
  
## INDEX FILE
When spingo is run, the reference database is converted internally to an efficient indexed structure. This process can take several minutes depending on the kmer size chosen and the size of the database. By default, this index is created each time the program is run. The index can be saved for reuse by using the --write-index (-w) option, e.g.:  
`SPINGODIR/spingo --write-index -d SPINGODIR/database/RDP_11.2.species.fa -i path/to/reads.fa > path/to/results.out`  
  
The index file is saved alongside the raw fasta database in the database directory. Once this is done, any use of the same database and kmer size will reuse the cached index file, saving several minutes per run.  
  
An additional program, spindex can be found alongside spingo. This program can be used to pre-generate these index files. See the section [SPINDEX](#SPINDEX) at the end of this document for more information.  
  
  
## SPINGO PARAMETERS
```
--help (-h)                 Display the help message

--version (-v)              Display the program version

--kmersize (-k) <int>       Kmer size. 
                            <int> is an integer in the range [1,15] (default=8).

--bootstrap (-b) <int>      Number of bootstrap samples to use.
                            <int> is a positive integer (default=10).

--subsample (-s) <int>      The proportion of kmers to use for each bootstrap subsample.
                            <int> is an integer > 0 (default is equal to kmersize).

--processors (-p) <int>     The number of processor threads to use.
                            <int> is an integer > 0 (default = 1).

--writeindex (-w)           A flag indicating that the index should be saved for future use.
                            See the INDEX FILE section above.

--ambiguous (-a)            A flag indicating that the output should contain an extra column
                            listing all species associated with an 'AMBIGUOUS' assignment.

--database (-d) <path>      Location of the fasta format database.

--input (-i) <path>         Location of the input file.
```
  
## OUTPUT FORMAT
Spingo writes output in a plain text format with tab separated columns. The column desctiptions are as follows: 
```
Column 1:   QUERY   Query sequence identifier.
Column 2:   SCORE   Similarity score (i.e. the fraction of k-mers in the query sequence that are shared with the reference sequence).
Column 3:   L1      Taxonomic assignment at level 1 (e.g. clostridium group).
Column 4:   L1_BS   Bootstrap score of the level 1 assignment.
Column 5:   L2      Taxonomic assignment at level 2 (e.g. genus).
Column 6:   L2_BS   Bootstrap score of the level 2 assignment.
Column 7:   L3      Taxonomic assignment at level 3 (e.g. species).
Column 8:   L3_BS   Bootstrap score of the level 3 assignment.
Column 9:   LIST    (Optional) list of species responsible for an ambiguous assignment.
```
  
## RESULTS SUMMARY
An additional script, spingo_summary can be found alonside spingo and spindex. This can be used to create a convenient summary of the results from a spingo run. e.g.  
`spingo_summary path/to/results.out`  
  
Usage:  
`spingo_summary [-h] RESULTS_FILE [--level <1,2>] [--similarity <float>] [--threshold <float>] [--percent]`
  
Parameters:  
```
--help (-h)                 Display the help message

RESULTS_FILE                Path to the input file (the output from spingo)

--level (-l) <int>          Taxonomic level to summarize.
                            1 = Clostridium group
                            2 = Genus
                            3 = Species

--similarity (-s) <float>   Similarity score cutoff, in the range [0,1]. Default is 0.5.

--threshold (-t) <float>    Bootstrap cutoff for the current level [0,1]. Default is 0.8.

--percent (-p)              A flag indicating that the summary should be expressed as percentages instead of raw read counts.
```

## SPINDEX
An additional program, spindex is supplied which can be used to pre-create index files before running spingo. Invoking spindex with the --help parameter will display a short list of available options.   
An example command to create an index for the RDP_11.2.species.fa database with a kmer size of 8 and using 4 processor threads:  
`SPINGODIR/spindex -k 8 -p 4 -d SPINGODIR/database/RDP_11.2.species.fa`  
  
  
## ALTERNATIVE TAXONOMY
It is possible to edit the existing taxonomy.map (or make a new one) to create your own taxonomic labels. The format of this file is self explanatory.  
  
  
## 64bit or 32bit
By default, `SPINGODIR/spingo` and `SPINGODIR/spindex` are symolic links to the 32bit precompiled binaries located in `SPINGODIR/dist/32bit/`. The 32bit binaries should run on both 32bit and 64bit operating systems.  
  
However, if you wish to use the precompiled 64bit versions instead, you can update the symbolic links to point at the binaries located in `SPINGODIR/dist/64bit/`  
e.g. from within SPINGODIR  
```
ln -f -s dist/64bit/spingo spingo
ln -f -s dist/32bit/spindex spindex
```
  

## COMPILING FROM SOURCE
SPINGO should compile on any POSIX compliant operating system with the following pre-requisites installed:  
- g++ (tested with g++ >= 4.1.2)  
- BOOST development libraries (>= 1.33), specifically:
  - Boost::program_options
  - Boost::serialization
  - Boost::thread
  
For RedHat derived distributions (eg CentOS) these can be installed using the distributions boost-devel package
  
For Debian derived distributions (eg Ubuntu) these can be installed individually using the packages
- libboost-program-options-dev
- libboost-serialization-dev
- libboost-thread-dev
  
For other distributions, search your package maintainer for the appropriate packages.
  
Compilation is simple. Change to the `SPINGODIR/source` directory and execute:  
`make`  
This will compile the executables and update the symbolic links in SPINGODIR to use the compiled versions.  
  
If you wish to revert to using the original pre-compiled binaries,  
`make clean`  
from within the `SPINGODIR/source` directory will restore the original symbolic links.  
  
  
## COMPILATION ERRORS
If you are using a recent version of the boost development libraries (specifically 1.5x), compilation of SPINGO may fail with the message:  
`undefined reference to symbol '_ZN5boost6system15system_categoryEv'`  
and / or:  
`error adding symbols: DSO missing from command line`  
  
In this case, compile SPINGO using the following command from within the `SPINGODIR/source` directory:  
`make boost-fix=1`  
