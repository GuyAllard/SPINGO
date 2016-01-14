#! /usr/bin/env python2.7

#===========================================================================
# Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
# File: v_ripper.py
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


import sys
import re
import argparse
import signal

from Bio              import SeqIO
from Bio.SeqRecord    import SeqRecord
from Bio.Seq          import Seq
from multiprocessing  import Pool

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("ambicode", "common/ambicode.py")
load_src("utility", "common/utility.py")
import ambicode
import utility

#------------------------------------------------------------------------------

def permute(primer):
    
    # allow for a single insertion
    inserts = ["%sN%s" % (primer[:x], primer[x:]) for x in range(1, len(primer))]
    # or a single substitution (internal)
    subs = ["%sN%s" % (primer[:x], primer[x+1:]) for x in range(1, len(primer)-1)]
    # or a single deletion
    dels = ["%s%s" % (primer[:x], primer[x+1:]) for x in range(len(primer))]
    
    return "(%s)" % "|".join([ambicode.expand_sequence(s) for s in [primer]+inserts+subs+dels])

#------------------------------------------------------------------------------

def load_primers(primer_file, allow_mismatch):
    """read in forward and reverse primers from a fasta formatted file stream"""
    
    forward_primer = ""
    reverse_primer = ""

    print >> sys.stderr, "loading primer sequences from %s" % primer_file.name

    for record in SeqIO.parse(primer_file, "fasta"):
        record_id = record.id.upper()
        record_desc = record.description.upper()
        
        if "FORWARD" in record_id or "FORWARD" in record_desc:
            forward_primer = str(record.seq)
            
        elif "REVERSE" in record_id or "REVERSE" in record_desc:
            reverse_primer = str(record.seq)
            reverse_primer_rc = str(record.seq.reverse_complement())

    if forward_primer == "" or reverse_primer == "":
        print >> sys.stderr, "Could not find primer sequences"
        return None

    print >> sys.stderr, ("forward primer: %s \t reverse primer: %s" % 
                          (forward_primer, reverse_primer))
    
    if allow_mismatch:
        return (permute(forward_primer), permute(reverse_primer), permute(reverse_primer_rc))
        
    return (ambicode.expand_sequence(forward_primer), 
            ambicode.expand_sequence(reverse_primer), 
            ambicode.expand_sequence(reverse_primer_rc))

#------------------------------------------------------------------------------

def trim_sequence(seq, seq_id, sizes, logger):
    if sizes:
        seq_len = len(seq)
        (min_size, max_size) = sizes
        
        if seq_len < min_size:
            #logger.log("%s : V region too short (%i) - discarding" % (seq_id, seq_len))
            return None
            
        # truncate long sequences
        if seq_len > max_size:
            #logger.log("%s : V region too long (%i) - discarding" % (seq_id, seq_len))
            return None
    
    return seq
    
#------------------------------------------------------------------------------

def quality_filter(seq, seq_id, badchars, logger):
    """discard any sequences containing ambiguous nucleotides"""
    
    if badchars > -1:
        if len(re.findall(r'[rykmswbdhvnx-]', seq, re.IGNORECASE)) > badchars:
            #logger.log("%s : dirty sequence - discarding" % seq_id)
            return None
    
    return seq
        
#------------------------------------------------------------------------------

def process_record((record, primers, sizes, trim, badchars, require, logger)):
    """ Extracts a region of the sequence bounded by primer pairs, and
    filtered to a specific length.

    Returns a new record with the extracted sequence, or None if the sequence
    did not match.
    """
    # eliminate sequences which are poor quality
    seq = quality_filter(str(record.seq), record.id, badchars, logger)
    if not seq:
        return None
        
    # extract and trim the primer bounded region
    # look for and trim forward primer
    fp = re.search(primers[0], seq, re.IGNORECASE)
    if fp:
        if "f" in trim:
            seq = seq[fp.span()[1]:]
        else:
            seq = seq[fp.span()[0]:]
    elif "f" in require:
        #logger.log("%s : Forward primer not found in sequence - discarding" % record.id)
        return None
        
    # then the reverse primer
    rp = re.search(primers[1], seq, re.IGNORECASE)
    if rp:
        if "r" in trim:
            seq = seq[:rp.span()[0]]
        else:
            seq = seq[:rp.span()[1]]
    else:
        # try the reverse compliment
        rp = re.search(primers[2], seq, re.IGNORECASE)
        if rp:
            if "r" in trim:
                seq = seq[:rp.span()[0]]
            else:
                seq = seq[:rp.span()[1]]
        elif "r" in require:
            #logger.log("%s : Reverse primer not found in sequence - discarding" % record.id)
            return None
    
    # eliminate sequences which are too long or short
    seq = trim_sequence(seq, record.id, sizes, logger)
    if not seq:
        return None

    # sequence made the grade, hand it back
    return SeqRecord(Seq(seq), record.id, "", "")

#------------------------------------------------------------------------------

def process_file(input_file, num_processes, primers, sizes, trim, badchars, require, logger):
    """Generator yielding only records from 'input_file' which contain
    valid V-regions.
    V-regions are identified using a string search
    """
    
    pool = Pool(num_processes)
    records = ((rec, primers, sizes, trim, badchars, require, logger) for rec in SeqIO.parse(input_file, "fasta"))
    chunks = []
    done = False
    while not done:
        try:
            # add 1 line to chunks
            chunks += [records.next()]
        except StopIteration: # no more lines left to read
            done = True
            
        # send the chunks off to be processed
        if len(chunks) == num_processes or len(chunks) and done:
            results = pool.map(process_record, chunks)
            chunks = []
            for res in results:
                if res:
                    yield res
        
#------------------------------------------------------------------------------

def parse_args():
    """Process command line arguments"""
    
    description = ("Extract variable region from sequences. ")

    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("input", 
                        help="Input file name, '-' for standard input", 
                        metavar="IN_FILE",
                        type=lambda x: utility.parse_file_arg(parser, x, "r"))
    
    parser.add_argument("output", 
                        help="Output file name, '-' for standard output",
                        metavar="OUT_FILE",
                        type=lambda x: utility.parse_file_arg(parser, x, "w"))
                        
    parser.add_argument("primerfile",
                        help="fasta file containing forward and reverse primers",
                        metavar="PRIMER_FILE",
                        type=lambda x: utility.parse_file_arg(parser, x, "r", False))
                        
    parser.add_argument("-l", "--logfile", dest="logger", metavar="LOGFILE",
                        help="For logging changes",
                        default="v_ripper_%s.log" % utility.time_stamp(), 
                        type=lambda x: utility.simple_logger(x))
                        
    parser.add_argument("--trim", "-t",  default="",
                        choices=["fr","f","r"],
                        help=("""which conserved regions to trim.
                                'fr' = forward and reverse, 
                                'f' = forward only, 'r' = reverse only"""))
                                
    parser.add_argument("--require", "-r", default="",
                        choices=["fr","f","r"],
                        help=("""which conserved regions are required for the
                                 sequence to be retained.
                                'fr' = forward and reverse, 
                                'f' = forward only, 'r' = reverse only"""))
                              
    parser.add_argument("--sizes", "-s", default=None,
                        type=int, metavar="N", nargs=2,
                        help="""min and max lengths for the final sequence. 
                        If omitted, sequences will not be filtered by size.""")
                        
    parser.add_argument("--mismatch", "-m",
                        help="""Allow for a single mismatch (substitution/indel)
                             to the primer sequence""",
                        action="store_true")
                        
    parser.add_argument("--badchars", "-b",
                        help="""maximum allowable number of ambiguous characters
                                (-1 for no filtering of bad characters)""",
                        type = int, default=-1)
                        
    parser.add_argument("--numprocs", "-n",
                        help="""Number of processors to use""",
                        type=int, default=1, choices=xrange(1,16))
                        
    args = parser.parse_args()
    args.sizes = sorted(args.sizes) if args.sizes else None
    return args
    
#------------------------------------------------------------------------------

@utility.timeit
def main():
    """Make it so"""
    
    # process command line args
    args = parse_args()
    
    # get primer sequences
    primers = load_primers(args.primerfile, args.mismatch)
    if not primers:
        sys.exit()
        
    # extract
    SeqIO.write( 
        process_file( args.input, 
                      args.numprocs,
                      primers, 
                      args.sizes,
                      args.trim, 
                      args.badchars,
                      args.require,
                      args.logger ), args.output, "fasta" )
    

#------------------------------------------------------------------------------

if __name__ == '__main__':
    signal.signal(signal.SIGINT, lambda s,f:sys.exit(0))
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("v_ripper.py - interrupted by user")
