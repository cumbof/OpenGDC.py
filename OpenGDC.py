#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Oct 10, 2020'

import sys, os, time, yaml, utils
import argparse as ap
from pathlib import Path

def read_params():
    p = ap.ArgumentParser( description = ( 'The OpenGDC.py script to extract, extend, and convert public genomic '
                                            'data and clinical metadata from the Genomic Data Commons portal' ),
                           formatter_class = ap.ArgumentDefaultsHelpFormatter )
    p.add_argument( '--tumor', 
                    type = str,
                    help = 'Case-sensitive GDC Tumor (e.g. TCGA-BRCA)' )
    p.add_argument( '--datatype', 
                    type = str,
                    help = 'Case-sensitive Experimental Data Type (e.g. "Methylation Beta Value")' )
    p.add_argument( '--after', 
                    type = str,
                    help = 'Date time used to filter data that must be retrieved' )
    p.add_argument( '--download',
                    action = 'store_true',
                    default = False,
                    help = 'Download genomic data and/or clinical metadata' )
    p.add_argument( '--download_dir',
                    type = str,
                    help = 'Path to the folder in which data will be located after download' )
    p.add_argument( '--convert',
                    action = 'store_true',
                    default = False,
                    help = 'Convert genomic data and/or clinical metadata' )
    p.add_argument( '--convert_dir',
                    type = str,
                    help = 'Path to the folder in which the converted data will be located' )
    p.add_argument( '--matrix',
                    type = str,
                    help = 'Export converted files to a data matrix (it works for a limited set of data types only)' )
    p.add_argument( '--settings',
                    type = str,
                    default = './settings.yaml',
                    help = 'Path to the settings.yaml file' )
    p.add_argument( '--verbose',
                    action = 'store_true',
                    default = False,
                    help = 'Print messages to STDOUT' )
    p.add_argument( '-v', 
                    '--version', 
                    action = 'version',
                    version = 'decrypt.py version {} ({})'.format( __version__, __date__ ),
                    help = "Print the current decrypt.py version and exit" )
    return p.parse_args()

if __name__ == '__main__':
    t0 = time.time()
    # init params
    args = read_params()

    if args.verbose:
        print( "Loading settings" )
    # Load settings
    with open( args.settings ) as settings_stream:
        settings = yaml.load( settings_stream, Loader=yaml.FullLoader )
    # Fixing assets file paths with absolute paths
    settings[ "assets" ][ "gencode" ] = os.path.abspath( settings[ "assets" ][ "gencode" ] )
    settings[ "assets" ][ "ncbi" ][ "history" ] = os.path.abspath( settings[ "assets" ][ "ncbi" ][ "history" ] )
    settings[ "assets" ][ "ncbi" ][ "reference" ] = os.path.abspath( settings[ "assets" ][ "ncbi" ][ "reference" ] )
    settings[ "assets" ][ "hugo" ] = os.path.abspath( settings[ "assets" ][ "hugo" ] )

    # Init list of downloaded files
    downloaded = [ ]
    if args.download:
        if args.verbose:
            print( "Downloading {} data for {}".format( args.datatype, args.tumor ) )
            print( "Save to directory: {}".format( args.download_dir ) )
        # Create download directory if it does not exist
        if not os.path.exists( args.download_dir ):
            os.mkdir( args.download_dir )
        # Start downloading data
        downloaded = utils.download( args.tumor.upper(), args.datatype, args.download_dir, 
                                     after_datetime=args.after, settings=settings, verbose=args.verbose )
    else:
        # If the download is not enabled, search for files into the download directory
        if os.path.exists( args.download_dir ):
            downloaded = list( Path( args.download_dir ).glob( '*.*' ) )
    
    # If the download directory is empty, exit
    if not downloaded:
        if args.verbose:
            print( "No files found for conversion" )
        sys.exit( 1 )

    # Init list of converted files
    converted_filepaths = [ ]
    if args.convert:
        # Create convert directory if it does not exist
        if not os.path.exists( args.convert_dir ):
            os.mkdir( args.convert_dir )
        
        if args.verbose:
            print( "Files that will be converted: {}".format( len( downloaded ) ) )
            print( "Loading external assets" )
        # Load external resources if required
        resources = utils.load_resources( args.datatype, settings, verbose=args.verbose )
        
        if args.verbose:
            print( "Defining header schema" )
        # Write header.schema file with info about bed file columns
        utils.dump_schema( args.datatype, args.convert_dir )

        # Start converting files in downloaded list
        for filepath in downloaded:
            print( "Converting {}".format( filepath ) )
            # Iteratively update resources if required
            converted, outfilepath, resources = utils.convert( args.datatype, filepath, args.convert_dir, 
                                                  settings, resources=resources, verbose=args.verbose )
            if converted:
                converted_filepaths.append( outfilepath )
    else:
        # If the conversion is not enabled, search for files into the convert directory
        if os.path.exists( args.convert_dir ):
            converted_filepaths = list( Path( args.convert_dir ).glob( '*.bed' ) )
    
    if converted_filepaths:
        if args.matrix:
            # TODO feature under development
            if args.verbose:
                print( 'Feature under development' )
    
    # Print total elapsed time and exit
    t1 = time.time()
    print( 'Total elapsed time {}s\n'.format( int( t1 - t0 ) ) )
