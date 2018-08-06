#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``main.py`` **module description**:

This module is the main module that 
- parses the input arguments,
- launches the comparison for all pairs of source CDS and target gene,
- computes ortho groups
- writes the results in output files

command line: 
when you don't have sturucture file:
python 'pathto main.py'  -it file -sf 'path to cds fasta file'  -tf 'path to gene fasta file' -s2tf 'path to assiciation cds vs. gene file'  -op 'output path' -of list -m SFA -f 'Yes you want to SFA_G or No if not ' -c 'blast or splign where you don't have sturucture file'

when you have sturucture file: 
python 'pathto main.py'  -it file -sf 'path to cds fasta file'  -tf 'path to gene fasta file' -s2tf 'path to assiciation cds vs. gene file'  -sef  'path to structure file' -op 'output path' -of list -m SFA -f 'Yes you want to SFA_G or No if not ' 


 

moduleauthor:: Aïda Ouangraoua , Safa Jammali and Jean-David Aguilar
Université de Sherbrooke, Canada
for running questions email us at Safa.Jammali@USherbrooke.ca

2017-2018

"""

# import external bilio

import os
import argparse
import glob

# import project-specific file
from getdata import *
from compare import *
from write import *



#####################
### Main ############


def build_arg_parser():
    """
    This function parses and controls all parameters given by the user in th command line

    Parameters
    ----------


    Returns
    -------
    parser
        list of arguments
    """

    parser = argparse.ArgumentParser(description="Transcript aligner")
    parser.add_argument('-c', '--choice', help="used where there 's no exon extremities (structre) information, it can be blast or splign")
    
    parser.add_argument('-f', '--force', help="SFA_G, it can be Yes or No")
    parser.add_argument('-m', '--method', help = " SFA (required)")
    parser.add_argument('-it', '--inputType', help="Input type : file, id or name (required)")
    parser.add_argument('-gn', '--geneNumber', help="Gene number: pairwise or multiple (required)")
    
    #required for it = file"
    parser.add_argument('-sf', '--sourceFile', help="Source file name (required if --inputType = file)")
    parser.add_argument('-tf', '--targetFile', help="Target file name (required if --inputType = file)")
    parser.add_argument('-s2tf', '--source2TargetFile', help="Source to target file name (required if --inputType = file)")
    parser.add_argument('-sef', '--sourceExonFile', help="Source Exon file name (required if --inputType = file)")
        
    #required for it = id or name, and gn = pairwise"
    parser.add_argument('-s1', '--firstSpecies', help="First species common name (required if --inputType = file or name, and --geneNumber = pairwise)") #Example : human
    parser.add_argument('-s2', '--secondSpecies', help="Second species common name (required if --inputType = file or name, and --geneNumber = pairwise)") #Example : mouse
    
    
    #required for it = id, and gn = pairwise"
    parser.add_argument('-gid1', '--firstGeneId', help="First gene Ensembl Id (required if --inputType = id, and --geneNumber = pairwise)") #Exemple : ENSG00000105695
    parser.add_argument('-gid2', '--secondGeneId', help="Second gene Ensembl Id (required if --inputType = id, and --geneNumber = pairwise)") #Exemple : ENSMUSG00000036634
    
    #required for it = id, and gn = multiple"
    parser.add_argument('-gidlf', '--geneIdListFile', help="Gene Id list file name (required if --inputType = id, and --geneNumber = multiple)")
    
    #required for it = name"
    parser.add_argument('-g', '--gene', help="Gene common name (required if --inputType = name)") #Example : MAG
    
    #required for it = name, and gn = multiple"
    parser.add_argument('-slf', '--speciesListFile', help="Species list file name (required if --inputType = name, and --genenumber = multiple)")
    
    parser.add_argument('-op', '--outputPrefix', help="Output prefix (required)")

    parser.add_argument('-of', '--outputFormat', help="Output format : list or aln (required)")

    return parser

if __name__ == '__main__':
    """
    This is the main function, it takes all arguments given and launches
        all module of the project
        
    Parameters
    ----------


    Returns
    -------
    
    """
    files1=glob.glob('./sequences/genes/*')
    for file1 in files1:
        os.remove(file1)
    files2=glob.glob('./sequences/cds/*')
    for file2 in files2:
        os.remove(file2)
    files3 = glob.glob('./results/splign_results/*')
    for file3 in files3:
        os.remove(file3)
    files4 = glob.glob('./results/ident_results/*')
    for file4 in files4:
        os.remove(file4)
    files5 = glob.glob('./results/blast_Results/*')
    for file5 in files5:
        os.remove(file5)
    # parses the input arguments
    parser = build_arg_parser()
    args = parser.parse_args()
    method = args.method
    outputPrefix = args.outputPrefix
    outputFormat = args.outputFormat
    choice = args.choice
    if (method == None):
        print "Argument -m <method> is required"
    
    if (outputPrefix == None):
        print "Argument -op <outputprefix> is required"

    if (outputFormat == None):
        print "Argument -of <outputformat> is required : list or aln"
    if (choice != 'splign'):
        print "blast is used to extract structure information"

    if(outputPrefix != None and outputFormat != None):
        print "Retrieving input data..."
        sourceData,targetData = get_data(args)
        nbinitialSource = len(sourceData)
	if args.sourceExonFile:
                 
        	write_input_files(outputPrefix,sourceData,targetData)

        # first module launches the pairwise spliced comparison for all pairs of source CDS and target gene by given source CDS and target genemethod
        print "Comparing sequences..."
        
        
        comparisonResults = spliceAlignment(sourceData, targetData,  method , args.force, choice, outputPrefix)
     
   
        print ' end comparing'
        # computes ortho groups
        print "Computing orthology groups..."
        orthologyGroups = computeOrthology(sourceData,targetData,comparisonResults)
        print "Completing orthology groups..."
        extendedSourceData,extendedOrthologyGroups = completeOrthology(sourceData,targetData,comparisonResults,orthologyGroups)
        
        # write the results in output files
        print "Writting output files"
        
        mblocklist = 0
        writeOutfile(outputPrefix,outputFormat,extendedSourceData,targetData,
                     extendedOrthologyGroups,
                     comparisonResults,nbinitialSource,mblocklist)
        
        
        

