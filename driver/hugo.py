#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Oct 10, 2020'

import os, bz2

# Load the HUGO DB in memory
def get_symbol_entrez_map( filepath ):
    symbol2entrez = { }
    # Read the flat file line by line
    with bz2.open( filepath, 'rt' ) as hugo:
        next( hugo ) # Skip the first row
        for line in hugo:
            if line.strip():
                line_split = line.split( "\t" )
                # Retrieve the mapping between gene symbol and entrez id only
                symbol = line_split[ 1 ].strip().lower()
                entrez = line_split[ 18 ]
                symbol2entrez[ symbol ] = entrez
    return symbol2entrez

# Check whether a gene symbol is in HUGO DB
def get_entrez_from_symbol( symbol2entrez, symbol_from_gencode ):
    if symbol2entrez:
        if symbol_from_gencode.strip().lower() in symbol2entrez:
            return symbol2entrez[ symbol_from_gencode.strip().lower() ]
    return None
