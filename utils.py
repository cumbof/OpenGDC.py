#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Oct 10, 2020'

import os, requests

# Import supported GDC data parsers
import parser.methylation as methylation

# Load drivers for external assets
import driver.gencode as gencode
import driver.ncbi as ncbi
import driver.hgnc as hgnc

# Define the set of available parsers
PARSERS = {
    'Methylation Beta Value': methylation
}

# Download data from the Genomic Data Commons
# Try submitting the same request "repeat" times in case of bad response status
def retrieve( url, locate, params, repeat=0 ):
    recursion_count = 0
    while recursion_count <= repeat:
        try:
            response = requests.get( url, headers={"Content-Type": "application/json"} )
            response.raise_for_status()
            with open( locate, 'wb' ) as file:
                file.write( response.content )
            break
        except:
            # Remove empty file in case of bad response status
            if os.path.exists( locate ):
                os.unlink( locate )
            # Repeat
            recursion_count += 1

# Make a query to the Genomic Data Commons and format the response as JSON
# Try submitting the same request "repeat" times in case of bad response status
def query( url, params, repeat=0 ):
    # Make a query to the 'files' endpoint with payload
    recursion_count = 0
    query_response = { }
    while recursion_count <= repeat:
        try:
            response = requests.post( url, headers={"Content-Type": "application/json"}, json=params )
            response.raise_for_status()
            query_response = response.json()
            break
        except:
            # Repeat
            recursion_count += 1
    return query_response

# Download data from the Genomic Data Commons
# Make a query to the GDC "files" endpoint to retrieve the list of files available for a given tumor and data type
# For each of the hit, start downloading by calling the "retrieve" function on the GDC "data" endpoint
# Use "after_datetime" to select files created after a specified date
def download( tumor, datatype, download_dir, after_datetime=None, settings=None, verbose=False ):
    if settings is None:
        if verbose:
            print( "Missing settings" )
        return [ ]
    # Search for data
    # Prepare a payload
    # Define a list of attributes that must be retrieved
    fields = [
        "file_name",
        "file_id"
    ]
    fields = ",".join(fields)
    # Define a filter with multiple conditions
    # cases.project.project_id = tumor              # Get data related to a particular tumor only
    # files.data_type = datatype                    # Get files related to a particular experimental type
    # files.access = open                           # Get public accessible data only
    filters = {
        "op":"and",
        "content":[
            {
                "op":"=",
                "content":{
                    "field":"cases.project.project_id",
                    "value":tumor
                }
            },
            {
                "op":"=",
                "content":{
                    "field":"files.data_type",
                    "value":datatype
                }
            },
            {
                "op":"=",
                "content":{
                    "field":"files.access",
                    "value":"open"
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

    # Append filter on data creation datetime if after_datetime is specified
    if after_datetime:
        params[ "filters" ][ "content" ].append(
            {
                "op":">=",
                "content":{
                    "field":"files.created_datetime",
                    "value":[
                        after_datetime
                    ]
                }
            }
        )
    
    if verbose:
        print("Querying GDC")
    # Submit a query to GDC to retrieve the list of available data
    query_response = query( settings[ "gdc" ][ "searchurl" ], params,
                            repeat=settings[ "gdc" ][ "repeat" ] )
    if not query_response:
        # Raise error: unable to query GDC
        return [ ]
    
    # Extract info from query_response
    data_filepaths = [ ]
    if len( query_response[ "data" ][ "hits" ] ) > 0:
        for hit in query_response[ "data" ][ "hits" ]:
            file_uuid = hit[ "file_id" ]    # Get the file uuid
            file_name = hit[ "file_name" ]  # Get the file name
            # Download uuid by querying the 'data' endpoint
            data_url = '{}{}?related_files=true'.format( settings[ "gdc" ][ "downloadurl" ], file_uuid )
            # Save file as <file_uuid>_<file_name>
            # Append the file uuid in front of the file name to retrieve the aliquot uuid during the conversion process
            data_path = os.path.join( download_dir, '{}_{}'.format( file_uuid, file_name ) )
            if not os.path.exists( data_path ):
                if verbose:
                    print( "\tDownloading {}_{}".format( file_uuid, file_name ) )
                # Start retrieving data
                retrieve( data_url, data_path, params, repeat=settings[ "gdc" ][ "repeat" ] )
            if os.path.exists( data_path ):
                data_filepaths.append( data_path )
    
    # Return the list of downloaded files
    return data_filepaths

# Convert GDC data
def convert( datatype, filepath, convert_dir, settings, resources={ }, verbose=False ):
    if settings is None:
        if verbose:
            print( "Missing settings" )
        return False, None, resources
    
    # Run a specific parser according to the specified "datatype"
    if datatype in PARSERS:
        # It returns a boolean value as the process exit code and the converted file path
        # It also returns the resurce back in case of updated
        return PARSERS[ datatype ].convert( datatype, filepath, convert_dir, settings, 
                                            resources=resources, verbose=verbose )
    
    if verbose:
        print( "Unsupported GDC Data Type" )
    return False, None, resources

# Load external resources
# Paths to the resource files are defined in settings.yaml
def load_resources( datatype, settings, verbose=False ):
    resources = { }
    # Define the external resources required to convert the Methylation Beta Value data
    if datatype == "Methylation Beta Value":
        # Load data from external assets
        if verbose:
            print( "\tLoading Gencode local DB" )
        # Do not load Gencode
        # Gencode must be partially loaded while converting
        resources[ "Gencode" ] = { } 
        if verbose:
            print( "\tLoading NCBI local DB" )
        # Load both NCBI reference and history files
        resources[ "NCBI" ] = {
            "DB": ncbi.get_symbol_entrez_map( settings[ 'assets' ][ 'ncbi' ][ 'reference' ] ),
            "Deprecated": ncbi.get_deprecated_symbol_entrez_map( settings[ 'assets' ][ 'ncbi' ][ 'history' ] )
        }
        if verbose:
            print( "\tLoading HGNC local DB" )
        # Load HGNC database
        resources[ "HGNC" ] = hgnc.get_symbol_entrez_map( settings[ 'assets' ][ 'hgnc' ] )
    return resources

# Dump the header.schema with the definition of the fields in the converted files
def dump_schema( datatype, convert_dir ):
    # Run a specific parser according to the specified "datatype"
    if datatype in PARSERS:
        PARSERS[ datatype ].dump_schema( convert_dir )
