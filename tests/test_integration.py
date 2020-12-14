#!/usr/local/bin/python
'''
End-to-end integration test.
'''

import os, subprocess, filecmp, shutil, shlex, tempfile, csv
import numpy as np

# msgsteiner is a fixture object
def test_forest_integration(msgsteiner):
    '''
    Forest integration test.  Runs Forest on the example data in the a549
    directory and compares the output files with stored reference files.
    '''
    # msgsteiner is parsed using conftest.py
    assert msgsteiner is not None, 'Please provide path to msgsteiner using --msgpath option'
    
    # Forest inputs
    curr_dir = os.path.dirname(__file__)
    forest_conf = os.path.join(curr_dir, '..', 'example', 'a549', 'tgfb_forest.cfg') #provided config file
    phos_weights = os.path.join(curr_dir, '..', 'example', 'a549', 'Tgfb_phos.txt')

    edge_file = os.path.join(curr_dir, '..', 'data', 'iref_mitab_miscore_2013_08_12_interactome.txt') #interactom

    # Create a tmp directory for output
    forest_out = tempfile.mkdtemp()

    #Arbitrary value
    seed = 2
    
    try:
        forest_path = os.path.join(curr_dir, '..', 'scripts', 'forest.py')
    
        fcmd='python %s --prize=%s --edge=%s --conf=%s  --outpath=%s --msgpath=%s --seed=%s'%(
            forest_path, phos_weights, edge_file, forest_conf, forest_out, msgsteiner, seed)
        subprocess.call(shlex.split(fcmd), shell=False)	

        match, mismatch, errors = filecmp.cmpfiles(forest_out, os.path.join(curr_dir, 'integration_test_standard'), 
                                   ['result_augmentedForest.sif', 'result_dummyForest.sif', 'result_edgeattributes.tsv',
                                    'result_info.txt', 'result_nodeattributes.tsv', 'result_optimalForest.sif'], 
                                   shallow=False)

        if len(match) != 6:
            print('Mismatching files: ', mismatch)
            print('Errors: ', errors)
            assert 0, 'Not all Forest output files match'

    except IOError:
        print('IO error')
    finally:
        shutil.rmtree(forest_out)

def test_garnet_integration():
    '''
    Garnet integration test.  Runs Garnet on the example FOXA1 data in the a549
    directory and compares the output files with stored reference files.
    '''
    # Garnet config file
    curr_dir = os.path.dirname(__file__)
    example_dir = os.path.join(curr_dir, '..', 'example', 'a549')
    orig_conf = os.path.join(example_dir, 'tgfb_foxa1_garnet.cfg')

    # Create a temporary copy of the config file that updates the relative
    # file paths, which are relative to the example/a549 directory
    conf_file = tempfile.NamedTemporaryFile(suffix='.cfg', delete=False)
    try:
        with open(orig_conf) as orig_conf_f:
            # Copy each line in the original config file and update
            # lines that contain the path to a file
            for line in orig_conf_f:
                # Matches File= or file=
                if 'ile=' in line:
                    tokens = line.split('=')
                    updated_file = os.path.normpath(os.path.join(example_dir, tokens[1]))
                    conf_file.write('%s=%s\n' % (tokens[0], updated_file))
                else:
                    conf_file.write(line)
    finally:
        conf_file.close()

    # Create a tmp directory for output
    # It must be a subdirectory of a known directory so that the relative
    # paths to the motif gifs are identical on different machines
    garnet_out = tempfile.mkdtemp(dir=curr_dir)

    try:
        garnet_path = os.path.join(curr_dir, '..', 'scripts', 'garnet.py')

        garnet_cmd='python %s --outdir=%s %s' % (garnet_path, garnet_out, conf_file.name)
        subprocess.call(shlex.split(garnet_cmd), shell=False)

        # Test all of the Garnet output files
        # Do not test whether Garnet produces extra files besides those listed here
        # Test events_to_genes_with_motifsregression_results.tsv as a special
        # case that does not use exact string matching for floating point values
        output_files = ['events_to_genes.fsa',
                        'events_to_genes.xls',
                        'events_to_genes_with_motifs.pkl',
                        'events_to_genes_with_motifs.tgm',
                        'events_to_genes_with_motifs.txt',
                        'events_to_genes_with_motifs_geneids.txt',
                        'events_to_genes_with_motifs_tfids.txt',
                        'events_to_genes_with_motifsregression_results.html',
                        'events_to_genes_with_motifsregression_results_FOREST_INPUT.tsv']
        match, mismatch, errors = filecmp.cmpfiles(garnet_out, os.path.join(curr_dir, 'integration_test_standard'),
                                   output_files, shallow=False)

        match_table, mismatch_table, errors_table = cmp_garnet_table(garnet_out, os.path.join(curr_dir, 'integration_test_standard'), 'events_to_genes_with_motifsregression_results.tsv')
        match.extend(match_table)
        mismatch.extend(mismatch_table)
        errors.extend(errors_table)

        # Add 1 because events_to_genes_with_motifsregression_results.tsv is
        # not in the output_files list
        if len(match) != (len(output_files)+1):
            print('Mismatching files: ', mismatch)
            print('Errors: ', errors)
            assert 0, 'Not all Garnet output files match'

    except IOError:
        print('IO error')
    finally:
        # Remove temporary config file here because delete=False above
        os.remove(conf_file.name)
        shutil.rmtree(garnet_out)

def cmp_garnet_table(test_dir, standard_dir, filename):
    '''
    A file comparison that mimics filecmp.cmp but is tailored for the Garnet
    output tsv file events_to_genes_with_motifsregression_results.tsv.  It
    does not rely on exact string matches but rather uses numpy.isclose
    for all numeric values.
    
    INPUT:
    test_dir - the test output directory
    standard_dir - the directory with the expected output file
    filename - the filename of the file that should exist in both directories
    
    OUTPUT:
    (match, mismatch, errors) - like filecmp.cmpfiles, lists of filenames
    that match, did not match, or produced errors except here the lists will
    contain exactly 0 or 1 item
    '''
    # Possible return values
    match = ([filename], [], [])    
    mismatch = ([], [filename], [])    
    error = ([], [], [filename])    

    try:
        # Load both tsv files
        with open(os.path.join(test_dir, filename)) as test_file:
            reader = csv.DictReader(test_file, delimiter = '\t')
            test_contents = list(reader)
        with open(os.path.join(standard_dir, filename)) as standard_file:
            reader = csv.DictReader(standard_file, delimiter = '\t')
            standard_contents = list(reader)
    except IOError:
        return error

    # Confirm the number of rows is the same
    if len(test_contents) != len(standard_contents):
        return mismatch

    # Confirm all rows are the same using string matching for the motif is
    # and approximate float matching for the slope, p-value, and q-value
    for row in range(len(test_contents)):
        test_row = test_contents[row]
        standard_row = standard_contents[row]

        if len(test_row) != len(standard_row):
            return mismatch
        if test_row['Motif'] != standard_row['Motif']:
            return mismatch
        for key in ['Slope', 'p-val', 'q-val']:
            if not np.isclose(float(test_row[key]), float(standard_row[key]), rtol=0, atol=1e-10, equal_nan=True):
                return mismatch

    return match
