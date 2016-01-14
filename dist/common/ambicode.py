#! /usr/bin/env python2.7

#===========================================================================
# Program: SPINGO - SPecies level IdentificatioN of metaGenOmic amplicons.
# File: ambicode.py
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


# lookup table for translating the ambiguity codes
ambitable = {'A':'A',       'a':'a',
             'T':'T',       't':'t',
             'G':'G',       'g':'g',
             'C':'C',       'c':'c',
             'R':'[AG]',    'r':'[ag]',
             'Y':'[CTU]',   'y':'[ctu]',
             'K':'[GTU]',   'k':'[gtu]',
             'M':'[AC]',    'm':'[ac]',
             'S':'[CG]',    's':'[cg]',
             'W':'[ATU]',   'w':'[atu]',
             'B':'[CGTU]',  'b':'[cgtu]',
             'D':'[AGTU]',  'd':'[agtu]',
             'H':'[ACTU]',  'h':'[actu]',
             'V':'[ACG]',   'v':'[acg]',
             'N':'[ACGTU]', 'n':'[acgtu]',
             'X':'X',       'x':'x',
             '-':'-',       '-':'-'}

#------------------------------------------------------------------------------

def expand_base( base ):
    """Returns a regex type string representation of base 's'
    with IUPAC ambiguity codes expanded.

    e.g.
    expand_base( 'N' ) => '[ACGTU]'
    """
    try:
        return ambitable[base]
    except KeyError:
        return base

#------------------------------------------------------------------------------

def expand_sequence( seq ):
    """Returns a regex type string representation of sequence 's'
    with IUPAC ambiguity codes expanded.

    e.g.
    expand_sequence( 'TAGCNTT' ) => 'TAGC[ACGTU]TT'
    """
    return "".join( expand_base(base) for base in seq )
