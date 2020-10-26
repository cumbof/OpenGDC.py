#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Oct 21, 2020'

import os, requests, utils, xmltodict

# Define supported input file extensions
def supported_ext( ):
    return [ "xml" ]

# header.schema definition
def dump_schema( convert_dir ):
    pass

# Define the conversion procedure for the Clinical and Biospecimen Supplements data type
def convert( datatype, filepath, outdir, settings, resources={ }, verbose=False ):
    datatype = None
    if "org_clinical." in os.path.basename( filepath ):
        datatype = "clinical"
    elif "org_biospecimen." in os.path.basename( filepath ):
        datatype = "biospecimen"
    else:
        return False, None, resources

    # File uuid is prepended to the file name and it is separated from the original file name by an underscore
    file_uuid = os.path.basename( filepath ).split( '_' )[ 0 ]
    if verbose:
        print( "\tProcessing {}".format( file_uuid ) )
    metadata = { }
    # Load XML to dict
    with open( filepath ) as xmlfile:
        metadict = xmltodict.parse( xmlfile.read(), dict_constructor=dict )
        if datatype == "clinical":
            clinical = { "clinical__{}".format( "__".join( [ str(k) for k in keypath[:-1] ] ) ): value for ( keypath, value ) in keypaths( metadict ) }
            # Search for patient_uuid
            patient_uuid = "NA"
            for key in clinical:
                if key.lower().endswith( "bcr_patient_uuid" ):
                    patient_uuid = clinical[ key ]
                    break
            # Use the "resources" channel to pass the results out
            resources[ patient_uuid ] = clinical
            return True, None, resources
        elif datatype == "biospecimen":
            biospecimen = { "__".join( [ str(k) for k in keypath ] ): value for ( keypath, value ) in keypaths( metadict )  }
            samples = { }
            for keypathstr in biospecimen:
                keypath = keypathstr.split( "__" )
                # Last element is always "#text"
                if keypath[ -2 ].lower().endswith( "bcr_aliquot_uuid" ):
                    aliquots = [ biospecimen[ keypathstr ] ]
                    if keypath[ 4 ] in samples:
                        aliquots = samples[ keypath[ 4 ] ] + aliquots
                    samples[ keypath[ 4 ] ] = aliquots            
            for sample_id in samples:
                sample_type_id_pathstr = "__".join( [ 'bio:tcga_bcr', 'bio:patient', 'bio:samples', 
                                                     'bio:sample', str(sample_id), 'bio:sample_type_id', '#text' ] )
                tissue_status = None
                if sample_type_id_pathstr in biospecimen:
                    # Retrieve sample type id
                    tissue_status = get_tissue_status( biospecimen[ sample_type_id_pathstr ] )

                for aliquot_id in range( len( samples[ str(sample_id) ] ) ):
                    data_map = { }
                    aliquot_path = [ 'bio:tcga_bcr', 'bio:patient', 'bio:samples', 'bio:sample', str(sample_id), 
                                     'bio:portions', 'bio:portion', 'bio:analytes', 'bio:analyte', str(aliquot_id), 
                                     'bio:aliquots', 'bio:aliquot' ]
                    for position in reversed( list( range( 1, len( aliquot_path ) + 1 ) ) ):
                        path = aliquot_path[ : position ]
                        for keypathstr in biospecimen:
                            keypath = keypathstr.split( "__" )
                            if ( len( set( keypath ).intersection( set( path ) ) ) == len( path ) and 
                                    len( keypath ) == len( path ) + 2 ):
                                data_map[ "biospecimen__{}".format( "__".join( keypath[:-1] ) ) ] = biospecimen[ keypathstr ]
                    
                    if tissue_status is None:
                        # If sample type id not in biospecimen
                        # Try to retrieve it from the aliquot barcode
                        aliquot_barcode_pathstr = '__'.join( aliquot_path + [ "bio:bcr_aliquot_barcode", "#text" ] )
                        sample_type_id = biospecimen[ aliquot_barcode_pathstr ].split( "-" )[ 3 ][:2]
                        tissue_status = get_tissue_status( sample_type_id )
                    data_map[ "biospecimen__tissue_status" ] = tissue_status
                    resources[ samples[ str(sample_id) ][ aliquot_id ].lower() ] = data_map

            return True, None, resources
    return False, None, resources

# Get mapping <key path, value>
def keypaths( nested ):
    for key, value in nested.items():
        if isinstance( value, dict ):
            for subkey, subvalue in keypaths( value ):
                yield [key] + subkey, subvalue
        elif isinstance( value, list ):
            for position in range( len( value ) ):
                tmpkey = [key] + [position]
                for subkey, subvalue in keypaths( value[ position ] ):
                    yield tmpkey + subkey, subvalue
        else:
            if key == "#text":
                yield [key], value

def get_tissue_status( sample_type_id ):
    sample_type_id = int( sample_type_id )
    if ( sample_type_id > 0 and sample_type_id < 10 ) or sample_type_id == 40:
        return "tumoral"
    elif sample_type_id > 9 and sample_type_id < 15:
        return "normal"
    elif sample_type_id == 20:
        return "control"
    else:
        return "undefined"

def build_metadata( outdir, clinical, biospecimen, verbose=False ):
    for aliquot_uuid in biospecimen:
        if verbose:
            print( "\tBuilding {}".format( aliquot_uuid ) )
        with open( os.path.join( outdir, "{}.meta".format( aliquot_uuid ) ), 'w+' ) as meta:
            for key in sorted( biospecimen[ aliquot_uuid ] ):
                meta.write( "{}\t{}\n".format( key, biospecimen[ aliquot_uuid ][ key ] ) )
            patient_uuid = biospecimen[ aliquot_uuid][ "biospecimen__{}".format( 
                                '__'.join( [ 'bio:tcga_bcr', 'bio:patient', 'shared:bcr_patient_uuid' ] ) ) ]
            if patient_uuid in clinical:
                for key in sorted( clinical[ patient_uuid ] ):
                    meta.write( "{}\t{}\n".format( key, clinical[ patient_uuid ][ key ] ) )
