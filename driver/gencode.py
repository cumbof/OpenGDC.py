#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Oct 10, 2020'

import os, bz2

# Load the Gencode DB partially
def get_gencode_info_fromfile( gencode_db, region_name, region_type, gencode_data={ } ):
    if not gencode_data:
        gencode_data = {
            'gene': { },
            'exon': { },
            'transcript': { },
            'utr': { },
            'cds': { },
            'start_codon': { },
            'stop_codon': { }
        }
    
    # If gencode_data is already defined, avoid reading the Gencode DB again
    if not gencode_data[ region_type.lower() ]:
        # Read the Gencode DB from assets
        with bz2.open( gencode_db, 'rt' ) as gencode:
            for line in gencode:
                line = line.strip()
                if line:
                    if not line.startswith( "#" ):
                        entry = { }
                        line_split = line.split( "\t" )
                        line_type = line_split[ 2 ].strip()
                        # Consider a specific type only
                        if line_type.lower() == region_type.lower():
                            # Retrieve chromosome, start, end, strand and type
                            entry[ 'chr' ] = line_split[ 0 ].strip()
                            entry[ 'start' ] = line_split[ 3 ].strip()
                            entry[ 'end' ] = line_split[ 4 ].strip()
                            entry[ 'strand' ] = line_split[ 6 ].strip()
                            entry[ 'type' ] = line_type
                            # Retrieve additional information from Gencode
                            extended_info_arr = line_split[ 8 ].strip().split( ";" )
                            symbol = "NA"
                            ensembl_id_noversion = "NA"
                            for data in extended_info_arr:
                                if data.lower().strip().startswith( "gene_name" ):
                                    symbol = data.strip().split( "\"" )[ -2 ]
                                elif data.lower().strip().startswith( "gene_id" ):
                                    ensembl_id_noversion = data.strip().split( "\"" )[ -2 ].split( "." )[ 0 ];
                            entry[ 'symbol' ] = symbol
                            entry[ 'ensembl_id' ] = ensembl_id_noversion
                            # Define an identifier for the Gencode map
                            identifier = ""
                            if region_name.lower() == "symbol":
                                identifier = symbol.lower()
                            elif region_name.lower() == "ensembl_id":
                                identifier = ensembl_id_noversion
                            # Put extended results into the Gencode map
                            if identifier.strip():
                                entries = [ ]
                                if identifier in gencode_data:
                                    entries = gencode_data[ region_type.lower() ][ identifier ]
                                entries.append( entry )
                                gencode_data[ region_type.lower() ][ identifier ] = entries
    return gencode_data