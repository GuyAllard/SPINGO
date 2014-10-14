#! /usr/bin/env python

#===========================================================================
# Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
# File: make_database.py
# Author: Guy Allard
# Copyright (C) 2014  Department of Microbiology, University College Cork 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================

from __future__ import print_function
import os.path
import sys
import argparse

#------------------------------------------------------------------------------

def open_input_file( parser, file_name):
    """Try and open the specified input file"""
    if not os.path.exists(file_name):
        parser.error("File '%s' does not exist" % file_name)
    else:
        try:    
            return open(file_name, "r")
        except IOError:
            parser.error("Could not open file '%s'" % file_name)

#------------------------------------------------------------------------------

class DnaSequence:
    def __init__(self, header="", sequence=""):
        self.header = header
        self.sequence = sequence

#------------------------------------------------------------------------------

def fasta_iterator( fasta_file ):
    seq = None
    for line in fasta_file:
        if line[0] == '>':
            if seq:
                yield seq
            seq = DnaSequence(line.strip()[1:])
        else:
            seq.sequence += line.strip()
        
    if seq:
        yield seq

#------------------------------------------------------------------------------

def parse_args():
    """Process command line arguments"""
    description = ("""Create database from RDP data""")

    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("taxonomy",
                       help="Taxonomy file name",
                       metavar="TAX_FILE",
                       type=lambda x: open_input_file(parser, x))
    
    return parser.parse_args()
    
#------------------------------------------------------------------------------

def main():
    # handle command line params
    args = parse_args()
     
    # load the taxonomy dictionary
    taxonomy_map = dict()
    for line in args.taxonomy:
        line = [x.strip() for x in line.split("\t")]
        taxonomy_map[line[0]] = "\t".join(line[1:])

    # process the file
    for record in fasta_iterator(sys.stdin):
        seq_id = [x.strip() for x in record.header.split()][0]
        try:
            tax = taxonomy_map[seq_id]
        except KeyError:
            continue

        print(">%s\t%s\n%s" % (seq_id, tax, record.sequence.upper()), file=sys.stdout)
    
#------------------------------------------------------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("make_database.py - interrupted by user")
