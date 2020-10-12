#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Oct 10, 2020'

import os, bz2

# Load the list of gene symbols currently in use from NCBI
def get_symbol_entrez_map( reference_filepath ):
    symbol2entrez = { }
    # Read the NCBI history line by line
    with bz2.open( reference_filepath, 'rt' ) as ref:
        for line in ref:
            line = line.strip()
            if line:
                if not line.startswith( "#" ):
                    extended_info = line.split( "\t" )[ 8 ]
                    if "Name=" in extended_info:
                        extended_info_arr = extended_info.split( ";" )
                        # Get the couple <deprecated gene symbol, entrez id>
                        symbol = ""
                        entrez = ""
                        for data in extended_info_arr:
                            if data.strip().lower().startswith( "name" ):
                                symbol = data.split( "=" )[ -1 ]
                            elif "geneid" in data.strip().lower():
                                for part in data.split( "," ):
                                    entrez = part.split( ":" )[ -1 ]
                                    break
                        symbol2entrez[ symbol.strip().lower() ] = entrez.strip().lower()
    return symbol2entrez

# Load the list of deprecated gene symbols from NCBI
def get_deprecated_symbol_entrez_map( history_filepath ):
    deprecated_symbol2entrez = { }
    # Read the NCBI history line by line
    with bz2.open( history_filepath, 'rt' ) as history:
        next( history ) # Skip the first row
        for line in history:
            line = line.strip()
            if line:
                line_split = line.split( "\t" )
                taxa_id = line_split[ 0 ]
                # Select rows which refer to the taxa id 9606 (Homo Sapiens)
                if taxa_id == "9606":
                    # Get the couple <deprecated gene symbol, entrez id>
                    deprecated_symbol = line_split[ 3 ].strip().lower()
                    entrez = line_split[ 2 ]
                    deprecated_symbol2entrez[ deprecated_symbol ] = entrez
    return deprecated_symbol2entrez

# Check whether a gene symbol is in HUGO DB
def get_entrez_from_symbol( symbol2entrez, deprecated_symbol2entrez, symbol_from_gencode ):
    # Search in reference
    if symbol2entrez:
        if symbol_from_gencode.strip().lower() in symbol2entrez:
            return symbol2entrez[ symbol_from_gencode.strip().lower() ]
    # Search for deprecated symbols
    if deprecated_symbol2entrez:
        if symbol_from_gencode.strip().lower() in deprecated_symbol2entrez:
            return deprecated_symbol2entrez[ symbol_from_gencode.strip().lower() ]
    return None
