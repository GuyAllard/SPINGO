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

BACT=release11_2_Bacteria_unaligned.fa.gz
ARCH=release11_2_Archaea_unaligned.fa.gz
TAX=taxonomy.map
OUT=RDP_11.2.species.fa

all: database

database: $(BACT) $(ARCH) $(TAX)
	@echo Generating database file $(OUT)
	zcat $(BACT) $(ARCH) | ../dist/make_database.py $(TAX) > $(OUT)
	@echo Done

$(BACT):
	@read -p "$(BACT) is missing. Download it from from the Ribosomal Database Project? (y/n): " dl; \
    if [ $$dl = "y" ]; then \
        wget http://rdp.cme.msu.edu/download/$(BACT); \
    else \
        echo; echo Manually place the $(BACT) file in this directory and re-run make; echo; \
        exit 1; \
    fi

$(ARCH):
	@read -p "$(ARCH) is missing. Download it from the Ribosomal Database Project? (y/n): " dl; \
    if [ $$dl = "y" ] || [ $$dl == "Y" ]; then \
        wget http://rdp.cme.msu.edu/download/$(ARCH); \
    else \
        echo; echo Manually place the $(ARCH) file into this directory and re-run make; echo;\
        exit 1; \
    fi
    
$(TAX):
	@echo; echo MISSING TAXNOMY FILE!; echo
	exit 1

