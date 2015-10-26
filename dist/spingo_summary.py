#!/usr/bin/env python

#===========================================================================
# Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
# File: spingo_summary.py
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

from __future__ import print_function, division
from collections import defaultdict
import os.path
import sys
import argparse

#-------------------------------------------------------------------------------

def open_input_file( parser, file_name):
    """Try and open the specified input file"""
    if not os.path.exists(file_name):
        parser.error("File '%s' does not exist" % file_name)
    else:
        try:
            return open(file_name, "r")
        except IOError:
            parser.error("Could not open file '%s'" % file_name)

#-------------------------------------------------------------------------------

def percentFloat ( string ):
    """Check that a given argument is within the range [0,1]"""
    value = float( string )
    if value < 0 or value > 1:
        raise argparse.ArgumentTypeError( 'Value has to be between 0 and 1' )
    return value

#-------------------------------------------------------------------------------

def parse_args():
    """Process command line arguments"""
    description = ("""Summarize spingo output""")

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("input", 
                        help="Input file name", 
                        metavar="RESULTS_FILE",
                        type=lambda x: open_input_file(parser, x))

    parser.add_argument("--level", "-l", default=3,
                        help="Level to summarize. Default=3 (species)",
                        type=int)

    parser.add_argument("--similarity", "-s",
                        help="Similarity score threshold [0-1]. Default=0.5",
                        type=percentFloat, default=0.5, metavar="N")

    parser.add_argument("--threshold", "-t",
                        help="bootstrap threshold cutoff [0-1]. Defualt=0.8",
                        type=percentFloat, default=0.8, metavar="N")

    parser.add_argument("--percent", "-p",
                        help="Display summary as a percentage instead of raw counts.",
                        action="store_true")

    return parser.parse_args()
 
#-------------------------------------------------------------------------------

def main():
    """Make it so"""
    
    # process command line args
    args = parse_args()
    
    if args.level <= 0:
        sys.exit("Error: level must be greater than 0")
    
    scoreField = 1
    assignmentField = args.level * 2
    bsField = args.level * 2 + 1
    
    results = defaultdict(int)
    total = 0

    # read in results
    for line in args.input:
        fields = line.strip().split("\t")
        if fields[0][0] == "#":
            continue

        if len(fields) < 2 + args.level * 2:
            sys.exit("Error: number of fields is less than the specified level")

        total += 1
        if (float(fields[bsField]) >= args.threshold and float(fields[scoreField]) >= args.similarity) or fields[assignmentField] == "AMBIGUOUS":
            results[fields[assignmentField]] += 1
        else:
            results["UNCLASSIFIED"] += 1    
    
    # output summary, highest to lowest, unclassified at end
    totalUnclassified = results["AMBIGUOUS"] + results["UNCLASSIFIED"]
    for key in sorted(results, key = lambda x: results[x], reverse=True):
        if key not in ["AMBIGUOUS", "UNCLASSIFIED"]:
            if args.percent:
                print ("%s\t%f" % (key, 100 * results[key] / total), file=sys.stdout)
            else:
                print ("%s\t%i" % (key, results[key]), file=sys.stdout)

    if args.percent:
        print ("UNCLASSIFIED\t%f\n(AMBIGUOUS\t%f)" % (100 * totalUnclassified / total, 100 * results["AMBIGUOUS"] / total), file=sys.stdout)
    else:
        print ("UNCLASSIFIED\t%i\n(AMBIGUOUS\t%i)" % (totalUnclassified, results["AMBIGUOUS"]), file=sys.stdout)
    
#-------------------------------------------------------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("spingo_summary.py - interrupted by user")
