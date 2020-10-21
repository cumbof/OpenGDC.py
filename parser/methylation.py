#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Oct 21, 2020'

import os, requests, utils
import driver.gencode as gencode
import driver.ncbi as ncbi
import driver.hgnc as hgnc

# Define supported input file extensions
def supported_ext( ):
    return [ "txt" ]

# header.schema definition
def dump_schema( convert_dir ):
    with open( os.path.join( convert_dir, 'header.schema' ), 'w+' ) as schema:
        schema.write(
            '<?xml version="1.0" encoding="UTF-8"?>\n'\
            '<gmqlSchemaCollection xmlns="http://genomic.elet.polimi.it/entities" name="GLOBAL_SCHEMAS">\n'\
                '\t<gmqlSchema type="tab" coordinate_system="1-based">\n'\
                    '\t\t<field type="STRING">chrom</field>\n'\
		            '\t\t<field type="LONG">start</field>\n'\
		            '\t\t<field type="LONG">end</field>\n'\
		            '\t\t<field type="CHAR">strand</field>\n'\
		            '\t\t<field type="STRING">composite_element_ref</field>\n'\
                    '\t\t<field type="DOUBLE">beta_value</field>\n'\
                    '\t\t<field type="STRING">gene_symbol</field>\n'\
                    '\t\t<field type="STRING">entrez_gene_id</field>\n'\
                    '\t\t<field type="STRING">gene_type</field>\n'\
                    '\t\t<field type="STRING">ensembl_transcript_id</field>\n'\
                    '\t\t<field type="STRING">position_to_tss</field>\n'\
                    '\t\t<field type="STRING">all_gene_symbols</field>\n'\
                    '\t\t<field type="STRING">all_entrez_gene_ids</field>\n'\
                    '\t\t<field type="STRING">all_gene_types</field>\n'\
                    '\t\t<field type="STRING">all_ensembl_transcript_ids</field>\n'\
                    '\t\t<field type="STRING">all_positions_to_tss</field>\n'\
                    '\t\t<field type="STRING">cgi_coordinate</field>\n'\
                    '\t\t<field type="STRING">feature_type</field>\n'\
                '\t</gmqlSchema>\n'\
            '</gmqlSchemaCollection>'
        )

# Define the conversion procedure for the Methylation Beta Value data type
def convert( datatype, filepath, outdir, settings, resources={ }, verbose=False ):
    # File uuid is prepended to the file name and it is separated from the original file name by an underscore
    file_uuid = os.path.basename( filepath ).split( '_' )[ 0 ]
    if verbose:
        print( "\tProcessing {}".format( file_uuid ) )
    # Prepare a payload
    # aliquot_id is the field that must be retrieved
    fields = "cases.samples.portions.analytes.aliquots.aliquot_id"
    # Define a filter with multiple conditions
    # files.file_id = file_uuid                 # Search for a specific file with the specified file_uuid
    # files.data_type = datatype                # Specify a data type
    filters = {
        "op":"and",
        "content":[
            {
                "op":"=",
                "content":{
                    "field":"files.file_id",
                    "value":file_uuid
                }
            },
            {
                "op":"=",
                "content":{
                    "field":"files.data_type",
                    "value":datatype
                }
            }
        ]
    }
    params = {
        "filters": filters,
        "fields": fields,
        "format": "JSON",                               # Set a response output format
        "size": str( settings[ "gdc" ][ "size" ] )      # Set the maximum amount of entries in the result
    }
    
    # Submit a query to GDC to retrieve the aliquot uuid
    query_response = utils.query( settings[ "gdc" ][ "searchurl" ], params,
                                  repeat=settings[ "gdc" ][ "repeat" ] )
    if not query_response:
        # Unable to retrieve aliquot_uuid
        return False, None, resources
    
    # Retrieve aliquot_uuid
    aliquot_uuid = query_response[ "data" ][ "hits" ][0][ "cases" ][0][ "samples" ][0][ "portions" ][0][ "analytes" ][0][ "aliquots" ][0][ "aliquot_id" ]
    
    # Take all the converted lines in memory
    # Then sort them by chromosome and genomic coordinates
    # Finally, dump them to the output file
    dataMapChr = { }

    # Open the input file
    with open( filepath ) as gdc:
        next( gdc ) # Skip header
        for line in gdc:
            line = line.strip()
            if line:
                # Fields are separated by tabs
                line_split = line.split( "\t" )
                chromosome = line_split[ 2 ]
                beta_value = line_split[ 1 ]
                gene_symbols_comp = line_split[ 5 ]
                if chromosome != "*" and beta_value.lower() != "na" and ( gene_symbols_comp.strip() != "" and gene_symbols_comp.strip() != "." ):
                    start = line_split[ 3 ]
                    end = line_split[ 4 ]
                    composite_element_ref = line_split[ 0 ]
                    gene_types_comp = line_split[ 6 ]
                    transcript_ids_comp = line_split[ 7 ]
                    positions_to_tss_comp = line_split[ 8 ]
                    cgi_coordinate = line_split[ 9 ]
                    feature_type = line_split[ 10 ]
                    
                    # Enxtend info by querying Gencode, NCBI, and HGNC
                    fieldsmap, resources = extract_fields( chromosome, gene_symbols_comp, start, end, gene_types_comp,
                                                           transcript_ids_comp, positions_to_tss_comp, settings, 
                                                           resources=resources )
                    strand = fieldsmap[ "strand" ]
                    gene_symbol = fieldsmap[ "symbol" ]
                    gene_type = fieldsmap[ "gene_type" ]
                    transcript_id = fieldsmap[ "transcript_id" ]
                    position_to_tss = fieldsmap[ "position_to_tss" ]
                    entrez_id = fieldsmap[ "entrez" ]
                    all_entrez_ids = fieldsmap[ "entrez_ids" ]
                    all_gene_symbols = fieldsmap[ "gene_symbols" ]
                    all_gene_types = fieldsmap[ "gene_types" ]
                    all_transcript_ids = fieldsmap[ "transcript_ids" ]
                    all_positions_to_tss = fieldsmap[ "positions_to_tss" ]

                    # Values in "values" list compose the output line
                    values = [ chromosome, start, end, strand, composite_element_ref, 
                               beta_value, gene_symbol, entrez_id, gene_type, transcript_id, 
                               position_to_tss, all_gene_symbols, all_entrez_ids, all_gene_types,
                               all_transcript_ids, all_positions_to_tss, cgi_coordinate, feature_type ]

                    # Build a nested dict with chromosome and start position as keys to sort lines by genomic coordinates
                    chromosome_id = int( chromosome.replace( "chr", "" ).replace( "X", "23" ).replace( "Y", "24" ) )
                    start_id = int( start )
                    dataMapStart = { start_id: values }
                    if chromosome_id in dataMapChr:
                        dataMapStart = dataMapChr[ chromosome_id ]
                        dataList = [ ]
                        if start_id in dataMapStart:
                            dataList = dataMapStart[ start_id ]
                        dataList.append( values )
                        dataMapStart[ start_id ] = dataList
                    dataMapChr[ chromosome_id ] = dataMapStart


    if dataMapChr:
        # The same aliquot can be used for multiple experiments
        # Add "-mbv" suffix to avoid conflicts
        bed_filepath = os.path.join( outdir, '{}-mbv.bed'.format( aliquot_uuid ) )
        with open( bed_filepath, 'w+' ) as bed:
            # Sort chromosomes
            chromosomes = sorted( dataMapChr.keys() )
            for chromosome in chromosomes:
                # Sort start positions
                start_positions = sorted( dataMapChr[ chromosome ].keys() )
                for start in start_positions:
                    # Finally dump the line out on the BED file
                    for entry in dataMapChr[ chromosome ][ start ]:
                        bed.write( '{}\n'.format( '\t'.join( [ str( value ) for value in entry ] ) ) )
        return True, bed_filepath, resources
    return False, None, resources

# Extract significant info and extend data by querying Gencode, NCBI, and HGNC
def extract_fields( chromosome, gene_symbols_comp, start_site, end_site, gene_types_comp,
                    transcript_ids_comp, positions_to_tss_comp, settings, resources={ } ):
    result = { }
    gene2CpGdistance = { }
    gene2DistanceFromCpG = { }
    gene2startEnd = { }
    
    # Define the list of info that must be retrieved
    transcript = ""
    position_to_TSS = ""
    gene_type = ""
    gene_symbol = ""
    strand = "*"
    entrez = ""
    all_entrez_ids = ""
    all_gene_symbols = ""
    all_gene_types = ""
    all_transcript_ids = ""
    all_positions_to_TSS = ""

    # From the OpenGDC Format Definition:
    #   gene_symbol (i.e., the symbol of each of the genes (can be more than one, separated by the ;
    #   char) associated with the CpG site. The association is defined with methylations whose region
    #   (2 bp) is superimposed (for at least 1 base) to the gene region (gene body) or to a neighborhood
    #   of 1,500 bp upstream of the gene. The same gene symbol is repeated if more than one
    #   transcript_id of the gene (reported in field 8) is associated with the methylation site.)
    # For each of the gene symbols defined in the original GDC data
    genes = gene_symbols_comp.split( ";" )
    for i in range( len( genes ) ):
        gene = genes[ i ]
        last = i + 1
        for pos in range( i + 1, len( genes ) ):
            if genes[ pos ] != genes[ i ]:
                last = pos
                break
        # Load gene info from Gencode
        resources[ "Gencode" ] = gencode.get_gencode_info_fromfile( settings[ "assets" ][ "gencode" ], 
                                                                    "symbol", "gene", gencode_data=resources[ "Gencode" ] )
        gencode_info = resources[ "Gencode" ][ 'gene' ][ gene.lower() ]
        if gencode_info:
            gene_info = gencode_info[ 0 ]
            if gene_info:
                # Retrieve start and end position of the current gene
                start = gene_info[ 'start' ]
                end = gene_info[ 'end' ]
                # Define the distance of the gene from the CpG site
                if int( start_site ) >= int( start ) and int( end ) >= int( end_site ):
                    distance = ( int( start_site ) - int( start ) ) + ( int( end ) - int( end_site ) )
                    gene2CpGdistance[ gene ] = distance
                else:
                    if int( end_site ) <= int( start ):
                        distance = int( start ) - int( end_site )
                    if int( end ) <= int( start_site ):
                        distance = int( start_site ) - int( end )
                    gene2DistanceFromCpG[ gene ] = distance
                
                # Retrieve the entrez gene id from NCBI
                # If not in reference, it could be deprecated
                entrez = ncbi.get_entrez_from_symbol( resources[ "NCBI" ][ "DB" ], 
                                                      resources[ "NCBI" ][ "Deprecated" ], 
                                                      gene )
                # If not in NCBI
                if entrez is None:
                    # Trying to retrive the entrez_id from GeneNames (HGNC) by gene symbol
                    entrez = hgnc.get_entrez_from_symbol( resources[ "HGNC" ], gene )
                    if entrez is None:
                        entrez = ""
        
        # Define the list of entrez ids, gene symbols, gene_types
        all_entrez_ids = '{};{}'.format( all_entrez_ids, entrez ) if all_entrez_ids.strip() else entrez
        all_gene_symbols = '{};{}'.format( all_gene_symbols, gene ) if all_gene_symbols.strip() else gene
        gene_type = gene_types_comp.split( ";" )[ i ]
        all_gene_types = '{};{}'.format( all_gene_types, gene_type ) if all_gene_types.strip() else gene_type
        index_entrez = [ i, last, entrez ]
        gene2startEnd[ gene ] = index_entrez
    
    # Retrieve the gene symbol which minimizes its distance from the CpG island
    if gene2CpGdistance:
        gene_symbol = min( gene2CpGdistance.items(), key=lambda x: x[ 1 ] )[ 0 ]
    else:
        if gene2DistanceFromCpG:
            gene_symbol = min( gene2DistanceFromCpG.items(), key=lambda x: x[ 1 ] )[ 0 ]
    
    if gene_symbol.strip():
        # Get start and end coordinates
        index_start = gene2startEnd[ gene_symbol ][ 0 ]
        index_end = gene2startEnd[ gene_symbol ][ 1 ]
        # Build the list of transcripts and positions to TSS
        for idx in range( index_start, index_end ):
            transcript = '{}|{}'.format( transcript, transcript_ids_comp.split( ";" )[ idx ] ) if transcript.strip() else transcript_ids_comp.split( ";" )[ idx ]
            position_to_TSS = '{}|{}'.format( position_to_TSS, positions_to_tss_comp.split( ";" )[ idx ] ) if position_to_TSS.strip() else positions_to_tss_comp.split( ";" )[ idx ]
        
        gene_type = gene_types_comp.split( ";" )[ index_start ]
        # Load gene info from Gencode 
        resources[ "Gencode" ] = gencode.get_gencode_info_fromfile( settings[ "assets" ][ "gencode" ], 
                                                                    "symbol", "gene", gencode_data=resources[ "Gencode" ] )
        gencode_info = resources[ "Gencode" ][ 'gene' ][ gene_symbol.lower() ]
        gene_info = gencode_info[ 0 ]
        # Get strand and entrez id
        strand = gene_info[ "strand" ]
        entrez = gene2startEnd[ gene_symbol ][ -1 ]

    # For each of the genes here involved
    for gene in gene2startEnd:
        index_start = gene2startEnd[ gene ][ 0 ]
        index_end = gene2startEnd[ gene ][ 1 ]
        transcript_tmp = ""
        position_to_TSS_tmp = ""
        for idx in range( index_start, index_end ):
            if transcript_tmp.strip():
                transcript_tmp += '|{}'.format( transcript_ids_comp.split( ";" )[ idx ] )
            else:
                transcript_tmp = transcript_ids_comp.split( ";" )[ idx ]
            if position_to_TSS_tmp.strip():
                position_to_TSS_tmp += '|{}'.format( positions_to_tss_comp.split( ";" )[ idx ] )
            else:
                position_to_TSS_tmp = positions_to_tss_comp.split( ";" )[ idx ]
        all_transcript_ids += ';{}'.format( transcript_tmp ) if all_transcript_ids.strip() else transcript_tmp
        all_positions_to_TSS += ';{}'.format( position_to_TSS_tmp ) if all_positions_to_TSS.strip() else position_to_TSS_tmp

    # Put extended fields in result
    result[ "gene_types" ] = all_gene_types
    result[ "gene_symbols" ] = all_gene_symbols
    result[ "entrez_ids" ] = all_entrez_ids
    result[ "transcript_ids" ] = all_transcript_ids
    result[ "positions_to_tss" ] = all_positions_to_TSS
    result[ "transcript_id" ] = transcript
    result[ "position_to_tss" ] = position_to_TSS
    result[ "gene_type" ] = gene_type
    result[ "strand" ] = strand
    result[ "symbol" ] = gene_symbol
    result[ "entrez" ] = entrez
    return result, resources
