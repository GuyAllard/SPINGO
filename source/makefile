#============================================================================\
Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.     \
File: makefile                                                               \
Author: Guy Allard                                                           \
Copyright (C) 2014  Department of Microbiology, University College Cork      \
                                                                             \
This program is free software: you can redistribute it and/or modify         \
it under the terms of the GNU General Public License as published by         \
the Free Software Foundation, either version 3 of the License, or            \
(at your option) any later version.                                          \
                                                                             \
This program is distributed in the hope that it will be useful,              \
but WITHOUT ANY WARRANTY; without even the implied warranty of               \
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                \
GNU General Public License for more details.                                 \
                                                                             \
You should have received a copy of the GNU General Public License            \
along with this program.  If not, see <http://www.gnu.org/licenses/>.        \
=============================================================================

CPP = g++
CPPFLAGS = -c -Wall
LFLAGS = -pthread -lboost_program_options -lboost_thread -lboost_serialization

ifdef boost-fix
	LFLAGS += -lboost_system
endif

all: release links

debug: CPPFLAGS += -g -DDEBUG
debug: spingo spindex

distrib: LFLAGS += -static
distrib: release

release: CPPFLAGS += -O3 -funroll-loops -DNDEBUG
release: spingo spindex

spingo: main.o fasta.o kmerizer.o database.o classifier.o
	$(CPP) main.o fasta.o kmerizer.o database.o classifier.o -o spingo $(LFLAGS)

spindex: makeindex.o fasta.o kmerizer.o database.o
	$(CPP) makeindex.o fasta.o kmerizer.o database.o -o spindex $(LFLAGS)

main.o: main.cpp
	$(CPP) $(CPPFLAGS) main.cpp

fasta.o: fasta.cpp
	$(CPP) $(CPPFLAGS) fasta.cpp

kmerizer.o: kmerizer.cpp
	$(CPP) $(CPPFLAGS) kmerizer.cpp

database.o: database.cpp
	$(CPP) $(CPPFLAGS) database.cpp

classifier.o: classifier.cpp
	$(CPP) $(CPPFLAGS) classifier.cpp

makeindex.o: makeindex.cpp
	$(CPP) $(CPPFLAGS) makeindex.cpp
	
clean:
	rm -rf *.o spingo spindex
	ln -f -s dist/32bit/spingo ../spingo
	ln -f -s dist/32bit/spindex ../spindex

links: spingo spindex
	ln -f -s source/spingo ../spingo
	ln -f -s source/spindex ../spindex
