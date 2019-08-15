# Written by Amanda Daigle
# Ernest Fraenkel's lab
# MIT Biological Engineering
# 2013-2014


import os
import sys
import argparse
from shutil import which

from OmicsIntegrator.forest import PCSFInput, PCSFOutput, crossValidation, \
    changeValuesAndMergeResults


def main():
    # Parsing arguments (run python PCSF.py -h to see all these decriptions)
    parser = argparse.ArgumentParser(
        description="Find multiple pathways within an interactome "
        "that are altered in a particular condition using the Prize Collecting"
        " Steiner Forest problem"
    )
    # required arguments
    parser.add_argument(
        "-p",
        "--prize",
        dest="prizeFile",
        help="(Required) Path to the text file "
        "containing the prizes. Should be a tab delimited file with lines:"
        '"ProteinName\tPrizeValue"',
    )
    parser.add_argument(
        "-e",
        "--edge",
        dest="edgeFile",
        help="(Required) Path to the text file "
        "containing the interactome edges. Should be a tab delimited file with"
        "3 or 4 columns: "
        '"ProteinA\tProteinB\tWeight(between 0 and 1)\tDirectionality(U '
        "or D, optional)",
    )
    # optional arguments
    parser.add_argument(
        "-c",
        "--conf",
        dest="confFile",
        help="Path to the text file containing "
        'the parameters. Should be several lines that looks like:'
        '"ParameterName ='
        ' ParameterValue". Must contain values for w, b, D.  May contain'
        ' values for optional parameters mu, garnetBeta, noise, r, g.'
        ' Default = "./conf.txt"',
        default="conf.txt",
    )
    parser.add_argument(
        "-d",
        "--dummyMode",
        dest="dummyMode",
        help="Tells the program which nodes "
        'in the interactome to connect the dummy node to. "terminals"= connect'
        ' to all terminals, "others"= connect to all nodes except for'
        ' terminals, "all"= connect to all nodes in the interactome. If'
        ' you wish you supply your own list of proteins, dummyMode could'
        ' also be the path to a text file containing a list of proteins'
        ' (one per line). Default = "terminals"',
        default="terminals",
    )
    parser.add_argument(
        "--garnet",
        dest="garnet",
        help="Path to the text file containing "
        "the output of the GARNET module regression. Should be a"
        '  tab delimited file with 2 columns: "TranscriptionFactorName'
        '\tScore". Default = "None"',
        default=None,
    )
    parser.add_argument(
        "--musquared",
        action="store_true",
        dest="musquared",
        help="Flag to add "
        "negative prizes to hub nodes proportional to their degree^2, rather"
        " than degree. Must specify a positive mu in conf file.",
        default=False,
    )
    parser.add_argument(
        "--excludeTerms",
        action="store_true",
        dest="excludeT",
        help="Flag to "
        "exclude terminals when calculating negative prizes. Use if you want"
        " terminals to keep exact assigned prize regardless of degree.",
        default=False,
    )
    parser.add_argument(
        "--outpath",
        dest="outputpath",
        help="Path to the directory which will "
        "hold the output files. Default = this directory",
        default=".",
    )
    parser.add_argument(
        "--outlabel",
        dest="outputlabel",
        help="A string to put at the beginning "
        'of the names of files output by the program. Default = "result"',
        default="result",
    )
    parser.add_argument(
        "--cyto30",
        action="store_true",
        dest="cyto30",
        help="Use this flag if "
        "you want the output files to be amenable with Cytoscape v3.0 (this is"
        " the default).",
        default=True,
    )
    parser.add_argument(
        "--cyto28",
        action="store_false",
        dest="cyto30",
        help="Use this flag if "
        "you want the output files to be amenable with Cytoscape v2.8, rather"
        " than v3.0.",
    )
    parser.add_argument(
        "--noisyEdges",
        dest="noiseNum",
        help="An integer specifying how many "
        "times you would like to add noise to the given edge values and re-run"
        " the algorithm. Results of these runs will be merged together and"
        ' written in files with the word "_noisyEdges_" added to their names.'
        'The noise level can be controlled using the configuration file.'
        ' Default = 0',
        type=int,
        default=0,
    )
    parser.add_argument(
        "--shuffledPrizes",
        dest="shuffleNum",
        help="An integer specifying how "
        "many times you would like to shuffle around the given prizes and"
        " re-run the algorithm. Results of these runs will be merged together"
        " and written in files with the word "
        '"_shuffledPrizes_" added to their names. Default = 0',
        type=int,
        default=0,
    )
    parser.add_argument(
        "--randomTerminals",
        dest="termNum",
        help="An integer specifying how many "
        "times you would like to apply your given prizes to random nodes in"
        "the interactome (with a similar degree distribution) and re-run the"
        " algorithm. Results of these runs will be merged together and written"
        ' in files with the word "_randomTerminals_" added to their'
        ' names. Default = 0',
        type=int,
        default=0,
    )
    parser.add_argument(
        "--knockout",
        dest="knockout",
        nargs="*",
        help="Protein(s) you "
        'would like to "knock out" of the interactome to simulate a knockout'
        ' experiment.',
        default=[],
    )
    parser.add_argument(
        "-k",
        "--cv",
        dest="cv",
        help="An integer specifying the k value if you "
        "would like to run k-fold cross validation on the prize proteins."
        " Default = None.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--cv-reps",
        dest="cv_reps",
        help="An integer specifying how many runs of "
        "cross-validation you would like to run. To use this option, you must"
        " also specify a -k or --cv parameter. Default = None.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "-s",
        "--seed",
        dest="seed",
        help="An integer seed for the pseudo-random "
        "number generators. If you want to reproduce exact results, supply"
        " the same seed. Default = None.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--merge",
        dest="merge",
        help="Don't merge results of multirun methods such as noisyEdges",
        default=False,
    )

    options = parser.parse_args()

    # Check cv parameters do not conflict
    if options.cv_reps is not None:
        if options.cv is None:
            sys.exit(
                "You cannot use the --cv-reps option without also specifying"
                " a k parameter for k-fold cross validation."
            )

    # Check if outputpath exists
    if not os.path.isdir(options.outputpath):
        sys.exit("Outpath %s is not a directory" % options.outputpath)

    # Ensure msgsteiner can be located before spending time parsing the input
    # files
    if which('msgsteiner') is None:
        sys.exit("ERROR: The msgsteiner code was not found on your path")
    # Process input, run msgsteiner, create output object, and write out
    # results
    inputObj = PCSFInput(
        options.prizeFile,
        options.edgeFile,
        options.confFile,
        options.dummyMode,
        options.knockout,
        options.garnet,
        options.shuffleNum,
        options.musquared,
        options.excludeT,
    )
    (edgeList, info) = inputObj.runPCSF(options.seed)
    outputObj = PCSFOutput(
        inputObj, edgeList, info, options.outputpath, options.outputlabel, 1
    )
    outputObj.writeCytoFiles(
        options.outputpath, options.outputlabel, options.cyto30
    )

    # Get merged results of adding noise to edge values
    if options.noiseNum > 0:
        merged = changeValuesAndMergeResults(
            'noisyEdges',
            options.seed,
            inputObj,
            options.noiseNum,
            options.outputpath,
            options.outputlabel,
            options.excludeT,
            merge=options.merge,
        )
        if merged is not None:
            merged.writeCytoFiles(
                options.outputpath,
                options.outputlabel + "_noisy",
                options.cyto30,
            )

    # Get merged results of shuffling prizes
    if options.shuffleNum > 0:
        merged = changeValuesAndMergeResults(
            'shufflePrizes',
            options.seed,
            inputObj,
            options.shuffleNum,
            options.outputpath,
            options.outputlabel,
            options.excludeT,
            merge=options.merge,
        )
        if merged is not None:
            merged.writeCytoFiles(
                options.outputpath,
                options.outputlabel + "_shuffled",
                options.cyto30,
            )

    # Get merged results of randomizing terminals
    if options.termNum > 0:
        merged = changeValuesAndMergeResults(
            'randomTerminals',
            options.seed,
            inputObj,
            options.termNum,
            options.outputpath,
            options.outputlabel,
            options.excludeT,
            merge=options.merge,
        )
        if merged is not None:
            merged.writeCytoFiles(
                options.outputpath,
                options.outputlabel + "_randomTerminals",
                options.cyto30,
            )

    # If k is supplied, run k-fold cross validation
    if options.cv is not None:
        if options.cv_reps is None:
            crossValidation(
                options.cv,
                1,
                inputObj,
                options.seed,
                options.outputpath,
                options.outputlabel,
            )
        else:
            for i in range(0, options.cv_reps):
                crossValidation(
                    options.cv,
                    i + 1,
                    inputObj,
                    options.seed,
                    options.outputpath,
                    options.outputlabel,
                )


if __name__ == "__main__":
    main()
