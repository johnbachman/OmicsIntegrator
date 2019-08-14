#!/usr/local/bin/python
'''
A test script to run through the sample dataset provided with garnet-forest
'''

__author__='Sara JG Gosline'
__email__='sgosline@mit.edu'

import os, sys
from optparse import OptionParser

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option('--forest-only',dest='forest_only',action='store_true',default=False,help='Set this flag to run forest on phospho-proteomic data only.  DEFAULT:%default')
    parser.add_option('--msgpath',dest='msgsteiner',type='string',help='Path to msgsteiner code, be sure to include!')
    parser.add_option('--doRandom',dest='rand',action='store_true',help='THIS WILL TAKE A LONG TIME: set this flag to do 50 permutations of forest using the --noisyEdges flag',default=False)
    
    opts,args=parser.parse_args()
    
    phos_weights='Tgfb_phos.txt'

    #garnet requires a configuration file that has all the data
    forest_out='tgfb_garnet_forest_output'
    garnet_conf='tgfb_garnet.cfg' #provided config file
    gcmd='python ../../scripts/garnet.py --outdir=%s %s'%(forest_out,garnet_conf) #command
    
    #forest requires more inputs
    forest_conf='tgfb_forest.cfg' #provided config file should include garnetBeta parameter
    edge_file='../../data/iref_mitab_miscore_2013_08_12_interactome.txt' #interactome

    msgsteinerpath=opts.msgsteiner ##WE NEED MSGSTEINER9 INSTALLED!!!


    
    ##now ready to run commands
    if not opts.forest_only:
        print(gcmd)
        res=os.system(gcmd)
        if res!=0:
            sys.exit('Error executing garnet, will not execute forest')
        garnet_output=forest_out+'/events_to_genes_with_motifsregression_results_FOREST_INPUT.tsv'
        #garnet_beta='0.1'
        fcmd='python ../../scripts/forest.py --prize=%s --edge=%s --conf=%s --garnet=%s --outpath=%s --msgpath=%s'%(phos_weights,edge_file,forest_conf,garnet_output,forest_out,msgsteinerpath)
        if opts.rand:
            fcmd=fcmd+' --noisyEdges=20'
        print('\n'+fcmd)
        os.system(fcmd)

    else:
        forest_out='tgfb_forest_output'
        if not os.path.exists(forest_out): ##FOREST WILL NOT CREATE DIRECTORY FOR YOU, GARNET WILL
            os.makedirs(forest_out)
        fcmd='python ../../scripts/forest.py --prize=%s --edge=%s --conf=%s  --outpath=%s --msgpath=%s'%(phos_weights,edge_file,forest_conf,forest_out,msgsteinerpath)
        if opts.rand:
            fcmd=fcmd+' --noisyEdges=50'
        print('\n'+fcmd)

        os.system(fcmd)

    
