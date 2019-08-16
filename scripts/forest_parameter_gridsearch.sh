#! /usr/bin/env bash
# Perform a grid search, running the steiner tree algorithm for different
# parameter combinations in parallel. Requires GNU parallel to be installed

programname=`basename "$0"`

function usage () {
    echo "usage: $programname [-p prizefile] [-e edgefile] [-o outpath]"
    echo "    -p    prizefile    prizefile for steiner run"
    echo "    -e    edgefile     edgefile for steiner run"
    echo "    -o    outpath      output directory to place results"
    echo "*********************************************************************"
    echo "Parameter lists for grid search and number of jobs can  "
    echo "be specified by editing this script. This script has GNU parallel as"
    echo "a dependency."
    echo "*********************************************************************"
    exit 1
}


###############################################################################
# Config. Edit for your own uses
# Number of parallel jobs to perform
njobs=1

# Parameter lists for grid search
w="1"
b="10"
D="10"
mu="0.001"
###############################################################################

while getopts 'p:e:o:h:' option
do
    case "$option" in
	p) p="$OPTARG";;
	e) e="$OPTARG";;
	o) o="$OPTARG";;
	h) usage; exit;;
    esac
done

if [ -z "$p"] || [ -z "$e"] || [ -z "$o"]
then
    usage
    exit 1
fi

outpath=$o/output_w{1}_b{2}_D{3}_mu{4}

parallel --no-notice -j $njobs run_forest.sh -w {1} -b {2} -d {3} -u {4} \
	 $p $e $outpath ::: $w ::: $b ::: $D ::: $mu






		
