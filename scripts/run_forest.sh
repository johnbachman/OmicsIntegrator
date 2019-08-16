#! /usr/bin/env bash

#########################################################
#Wrapper around the forest.py script in Omics Integrator
# Automates generation of temporary config files 
#Albert Steppi 7.27.2018
#! /usr/bin/env bash
#########################################################

programname=`basename "$0"`
DIR="$(dirname "$(readlink -f "$0")")"

function usage () {
    echo "usage: $programname [-w W] [-d Depth] [-b Beta] [-u Mu] [-r num_runs]"
    echo "                     [-n noise] [prize] [edge] [outpath]"
    echo "    -w    W        Max number of trees in generated forest"
    echo "    -d    Depth    Max path-length from dummy node to terminal nodes"
    echo "    -b    Beta     Controls tradeoff between including more terminals"
    echo "                   and using less reliable edges"
    echo "    -u    Mu       Controls downweighting of prizes based on degree."
    echo "                   Used to penalize hub nodes."
    echo "    -r    num_runs Number of noisy edge runs to generate. If unset,"
    echo "                   there will be no noisy edge runs"
    echo "    -n    noise    Standard deviation of Gaussian noise to add to"
    echo "                   edge weights for noisy runs."
    echo "          prize    Path to prize file."
    echo "          edge     Path to interactome."
    echo "          outpath  Directory to output results."
    echo "*********************************************************************"
    echo "Wrapper script around forest.py that generates config file based on"
    echo "command line arguments. Makes it easier to script calls to forest.py."
    echo "*********************************************************************"
    exit 1
}
# Default parameters
w=5 # controls the number of trees in output
d=10 # maximum depth from the dummy node
b=2 # controls the number of terminal nodes included
u=0 # penalize hubs with high degree
n=0.05 # std deviation of random noise to add to edge weights
# Optional r (number of noisy runs to perform)

# parse flagged arguments
while getopts 'w:d:b:u:r:n:' option
do
    case "$option" in 
	w) w="$OPTARG";;
	d) d="$OPTARG";;
	b) b="$OPTARG";;
	u) u="$OPTARG";;
	r) r="$OPTARG";;
	n) n="$OPTARG";;	
    esac
done

# positional arguments
prize=${@:$OPTIND:1}
edge=${@:$OPTIND+1:1}
outpath=${@:$OPTIND+2:1}

if [ -z "$prize" ] || [ -z "$edge" ] || [ -z "$outpath" ]
then
    echo $prize
    echo $edge
    echo $outpath
    usage
    exit 1
fi

# if the specified output directory does not exist, create it
mkdir -p $outpath

# Create a temporary folder for storing forest's config file
conf=$(mktemp -d "$TMPDIR/$(basename $0).XXXXXXXXXXXXX")
cat <<EOF>> $conf/forest_cfg
w = $w
D = $d
b = $b
mu = $u
noise = $n
EOF

echo $n
echo $([[ -n $n ]] && echo "--noisyEdges $n")
python $DIR/forest.py --prize=$prize --edge=$edge --outpath=$outpath --conf=$conf/forest_cfg $([[ -n $r ]] && echo "--noisyEdges $r")

# Add the generated conf file to the output folder so we can recover
# the parameters later
cp -r $conf/forest_cfg $outpath

# cleanup temporary files on exit
trap "{ rm -rf $conf; }" EXIT
