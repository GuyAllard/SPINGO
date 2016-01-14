#! /usr/bin/env python2.7

#===========================================================================
# Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
# File: utility.py
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
import gzip
import tarfile
import mimetypes
import argparse
import time

#------------------------------------------------------------------------------

def open_file( file_name, mode="r", allow_std=True ):
    """
    Opens a file. 
    Attempts to detect the file type in order to handle gzip and tar files.

    file_name -- the name of the file to open
    mode      -- the file access mode (eg "r" or "w")
    allow_std -- can this file be stdin/stdout
    """
    
    # don't do anything if file_name is not a string
    if not isinstance(file_name, basestring):
        return file_name
    
    # pick stdin/stdout if requested
    if file_name == "-":
        if allow_std:
            if mode == "r":
                return sys.stdin
            elif mode == "w":
                return sys.stdout
        else:
            raise IOError("stdin/stdout not allowed")
        
    # try and determine the file type
    file_type = mimetypes.guess_type(file_name)

    # tarfile
    if file_type[0] == "application/x-tar":
        return tarfile.open(file_name, mode)
    
    # gzip
    elif file_type[1] == "gzip":
        return gzip.open(file_name, mode)
    
    # normal file
    else:
        return open(file_name, mode)

#------------------------------------------------------------------------------

class simple_logger():
    """Simple message logger class"""
    def __init__(self, file_name = None):
        self.log_file = None
        if file_name is not None and file_name != "":
            try:
                self.log_file = open_file(file_name, "w", False)
            except IOError:
                print >> sys.stderr, "Could not open '%s' for logging" % file_name

    def __del__(self):
        if self.log_file:
            self.log_file.close()
            self.log_file = None
            
    def log(self, message):
        if self.log_file:
            print >> self.log_file, message

#------------------------------------------------------------------------------

def parse_file_arg(parser, file_name, mode, allow_std=True):
    """Tries to open a file from a command line argument"""
    try:
        return open_file(file_name, mode, allow_std)
    except IOError:
        parser.error("Could not open file '%s'" % file_name)
    
#------------------------------------------------------------------------------

def id_to_name(identifier):
    """returns the species name from the sequence identifier"""
    return "_".join(identifier.split("_")[1:3])
    
#------------------------------------------------------------------------------

def time_stamp():
    now = "%f" % time.time()
    now = now.split(".")
    now[1] = now[1][:3]
    return "".join(now)
    

#------------------------------------------------------------------------------

# decorator for checking execution time of functions
# from http://www.andreas-jung.com/contents/a-python-decorator-for-measuring-the-execution-time-of-methods
def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print '%r (%r, %r) %2.2f sec' % \
              (method.__name__, args, kw, te-ts)
        return result

    return timed
    
#------------------------------------------------------------------------------

def percentFloat ( string ):
    """Check that a given argument is within the range [0,1]"""
    value = float( string )
    if value < 0 or value > 1:
        raise argparse.ArgumentTypeError( 'Value has to be between 0 and 1' )
    return value
