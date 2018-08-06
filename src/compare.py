#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``compare.py`` **module description**:

This module is the module that launches the pairwise comparison for all pairs of source CDS and target gene
by given source CDS , target gene and the method of comparison.

moduleauthor:: Aïda Ouangraoua , Safa Jammali and Jean-David Aguilar
Université de Sherbrooke Canada

2017-2018

"""

from multiprocessing import Process, Queue
from multiprocessing import active_children

import os

import alignment
from alignment import *
import Exact
from  Exact import *
from write import *

# Constant Declaration
STATUS_EXISTING_PROTEIN = 1
STATUS_PREDICTED_PROTEIN = 2
STATUS_PREDICTED_CDS = 3
STATUS_PARTIAL_CDS = 5
STATUS_COMPLETE_CDS_DIFFERENT_STRUCTURE = 4
MIN_IDENTITY_FINAL = 0.5
MIN_ERROR = 50
GLOBAL_LIMIT = 10**4
CUTOFF_IDENTITY = 30
ALPHA = 21
EVALUE_STRUCTURE= 1.0/10000000 # 10-7
#########################
### COMPARISON ##########
#########################


def spliceAlignment(sourceData, targetData,  method,  force, choice, outputPrefix):

    """
    This function launch a parallel pairwise alignment between all CDS and all genes

    Parameters
    ----------
    sourceData: list
    	list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
    	list that contain list of all informations about all genes:  id and sequence of each gene

    method: method used to do comparison of CDS vs. gene.
           It can be SFA, Splign or Exact
   
    force: Yes or No
    	Yes to do the DP  to align whole remaining exons and No : not to do it
    choice: string
	splign or blast
    outputPrefix:string
	prefix for the output 

    Returns
    -------
    comparisonresults: list
        list of hits detected by the method of alignment
    """



    comparisonResults = []
    exonEndList = {}
    exonStartList = {}
    cdsExon = {}
    cds2GeneExon = {}
    geneExon = {}
    intronList = {}


    # put cds in fasta format file
    for cds in sourceData:
        cdsId, cdsSeq, null, null = cds
        cdsSeqFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
        cdsFile = open(cdsSeqFile, "w")
        cdsFile.write(">" + cdsId + "\n")
        cdsFile.write(cdsSeq)
        cdsFile.close()

    # put gene in fasta format file
    for gene in targetData:
        geneId, geneSeq = gene
        geneSeqFile = os.getcwd() + '/src/sequences/genes/Gene_' + geneId + '.fasta'
        geneFile = open(geneSeqFile, "w")
        geneFile.write(">" + geneId + "\n")
        geneFile.write(geneSeq)
        geneFile.close()


    for gene in targetData:
        geneId, geneSeq = gene
        exonStartList[geneId] = []
        exonEndList[geneId] = []
        geneExon[geneId] = []
        intronList[geneId] = []
    sourceData2=[]
    # Extract information about exons and introns for each gene based on theirs CDS
    for cds in sourceData:
	    blockList_ =[]
            cdsId, cdsSeq, cdsGeneId, blockList = cds
            cdsExon[cdsId] = []
            cdsSeqFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
            geneSeqFile = os.getcwd() + '/src/sequences/genes/Gene_' + cdsGeneId  + '.fasta'
	   
            if(len(blockList) == 0):
		
		if  choice == 'splign':
		        splignOutputName = launch_splign(cdsSeqFile, geneSeqFile, cdsId, cdsGeneId )
		        blockList_ = parse_splign_output(splignOutputName)
			blockList = blockList_
                       
		elif choice == 'blast':
			
			for gene in targetData:
        			geneId, geneSeq = gene
				if cdsGeneId == geneId:

					blockList_, blockList = BlastLocal(cdsSeqFile, cdsSeq, geneSeqFile,geneSeq,  cdsId, cdsGeneId, EVALUE_STRUCTURE )
					
					break
				else :
					pass
		else:
			print 'choose method to extractexon strusture by splign or blast'
           		exit(-1)
	    sourceData2.append([cdsId, cdsSeq, cdsGeneId, blockList_])
    

            intronList[cdsGeneId ] += extractIntronGene(cdsGeneId , cdsId, blockList)
            exonEndList[cdsGeneId ] += [block[3] for block in blockList]
            exonStartList[cdsGeneId ] += [block[2] for block in blockList]
            geneExon[cdsGeneId ] += [block[2:] for block in blockList]
           
            cdsExon[cdsId] = [block[:2] for block in blockList]

            verifexon(cdsId, cdsExon[cdsId])
            cds2GeneExon[cdsId] = [block[2:] for block in blockList]
   
    if len(sourceData2[0][3])!=0:
            
	    write_input_files(outputPrefix,sourceData2,targetData)
    else:
	pass
    
    # Order the blokcs list
    for gene in targetData:
        geneId, geneSeq = gene
        exonStartList[geneId].sort()
        exonEndList[geneId].sort()
        geneExon[geneId].sort()
        intronList[geneId].sort()


        # remove exons redundancy
        for i in range(len(geneExon[geneId]) - 1, 0, -1):
            if (geneExon[geneId][i] == geneExon[geneId][i - 1]):
                geneExon[geneId].remove(geneExon[geneId][i])
        # remove introns redundancy
        for i in range(len(intronList[geneId]) - 1, 0, -1):
            if (intronList[geneId][i] == intronList[geneId][i - 1]):
                intronList[geneId].remove(intronList[geneId][i])

    
    results = {}
    threads = []
    q = Queue()

    # launch parallel comparison for CDS against genes
    for gene in targetData:
        for cds in sourceData:

            t = Process(target=parallelCompare, args=(gene, cds, cdsExon, \
                                                      sourceData, targetData, exonStartList, exonEndList, cds2GeneExon,
                                                      intronList, geneExon, method,  force, q))

            threads.append(t)
            t.daemon = True
            t.start()

    # recovering results of all pairs
    while len(results.keys()) < len(targetData) * len(sourceData):
        item = q.get()
        if item is None:
            break
        else:
            results[item[5] + item[6]] = [item[0], item[1], item[2], item[3], item[4]]
   
    # delete CDS fasta files created before
    for cds in sourceData:
        cdsId, cdsSeq, null, null = cds
        cdssSeqFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
        os.remove(cdssSeqFile)
    
    # delete gene fasta files created before and put results in a list according to th first order of gene and CDS
    for gene in targetData:
        geneId, geneSeq = gene
        geneSeqFile = os.getcwd() + '/src/sequences/genes/Gene_' + geneId + '.fasta'
        os.remove(geneSeqFile)
        comparisonResults.append([])
        for cds in sourceData:
            cdsId, null, null, null = cds
            comparisonResults[-1].append(results[geneId + cdsId])
    
    # return list of blocks of alignment for each pair CDS vs. gene
    return comparisonResults


def verifexon(cdsId, SuccessivBlockList):
    '''
    This function help to know if cdsexon are successive (OK) or not (not valid)

    Parameters
    ---------- 
    cdsId:string
        cds Id
    SuccessivBlockList:list
    list of exon 

    Returns
    -------
    cut the program if the list of exons are not successive  

    '''
    for i in range (0,len(SuccessivBlockList)-1):
		
		if SuccessivBlockList[i][1]==  SuccessivBlockList[i+1][0]:
			pass
		else:
			print cdsId, '  cdsexon non valide , look: ', SuccessivBlockList, '\n'
			exit(-1)

def extractIntronGene(geneId, cdsId, blockList):
    """
    This function extracts introns list from gene based on its cds exons list given in blocklist.
    The blockList contains the list of intervals of the exon

    Parameters
    ----------

    geneId: string
        gene Id
    cdsId:string
        cds Id
    blockList: list
        list of blocks


    Returns
    -------
    intronList:list
        list of intron
    """


    exonList = ''
    coord = []
    intronList = []
    for block in blockList:
        exonList = block[2:]
        coord.append(exonList[0])
        coord.append(exonList[1])

    for i in range(1, len(coord) - 2, 2):
        beginIntron = coord[i]
        endIntron = coord[i + 1]
        intronList.append([beginIntron, endIntron])

    return intronList


def parallelCompare(gene, cds, cdsExon,  sourceData, targetData, exonStartList, exonEndList, cds2GeneExon,intronList,\
                    geneExon,method, force,  q):
    """
    This function launchs comparison according to the method choosed to do the alignment and returns analysed results

    Parameters
    ----------

    gene: list
            list contains gene Id and gene sequence

    cds: list
        list contains CDS Id, CDS sequence and the gene of this CDS
    cdsExon: dictionary
        dictionay with CDS id as key and exons interval list as value
    sourceData: list
        list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
        list that contain list of all informations about all genes:  id and sequence of each gene
    exonStartList: dictionary
        dictionay with gene id as key and exons start list as value
    exonEndList: dictionary
        dictionay with gene id as key and exons end list as value

    cds2GeneExon: dictionary
        dictionay with cds id as key and exons blocks list as values
    intronList:dictionary
               dictionary of gene intron, the keys are gene Id and the values are the blocks interval of intron gene
    geneExon: list
        list of exon of gene
    method: string
            mthod choosed to pairwise alignment

    force: Yes or No
    	Yes to do the DP  to align whole remaining exons and No : not to do it

    Returns
    -------
    q: list
        bloks and information about the alignment of each pairs of CDS against gene
    """



    cdsId, cdsSeq, cdsGeneId, null = cds
    cdsLen = len(cdsSeq)

    geneId, geneSeq = gene
    geneLen = len(geneSeq)

    print "Comparing CDS",cdsId,"with Gene", geneId
   
    if method == 'SFA':

        #tBlastx for local alignment
        blockListLocal  = localAlignment(cdsId,cdsLen, cdsExon[cdsId],geneId, geneLen, geneSeq)
      
        print 'finish SFA_L', cdsId, geneId, '\n'
      
        blockList=[]

        if  force == 'Yes':
		
	
		blockList = ExtendAndAlignRestOfExons(cdsId, geneId, blockListLocal, cdsLen, cdsSeq, cdsExon[cdsId], cdsExon, geneLen, geneSeq, geneExon, intronList)
		print 'finish SFA_G', cdsId, geneId, '\n'
        else:
		blockList = OnlyExtend(cdsId, geneId, blockListLocal, cdsLen, cdsSeq, cdsExon[cdsId], geneLen, geneSeq, geneExon, intronList)
		print 'finish SFA_E', cdsId, geneId, '\n'

    blockList = order(blockList)
    # analyse blocks results
    status, splicingSites, cdsTarget, texisting = analyse_result(blockList , cdsId, cdsSeq, geneId, geneSeq, sourceData, cds2GeneExon, cdsExon)
    # list to return
    q.put([status, blockList, splicingSites, cdsTarget, texisting, geneId, cdsId])

def OnlyExtend(cdsId, geneId, blockListLocal, cdsLen, cdsSeq, cdsExon, geneLen, geneSeq, geneExon, intronList):
    """
    This function extend the rest of exons segments that weren't aligned by the local alignement

    Parameters
    ----------

    geneId: string
  		gene Id
    cdsId: string
		cds Id
    blockListLocal: list
        list of blocks found by the local alignment
    cdsExon: list
		list of exons of cds
    geneLen: int
		gene length
    cdsLen: int
		cds length

    cdsSeq: string
	sequence of CDS 
    geneSeq: string
	sequence of gene

    geneExon:
    intronList:
   
    Returns
    -------
    blockList: list
        blocks from local alignment with blocks found by the extension

    """

    blockListextend = []
    blockListextended = []

    
    blockListextend = extractBlockListToextend(cdsId, geneId, blockListLocal, cdsExon, geneLen, cdsLen)
  
    if len(blockListextend) != 0:

        blockListextended = extendBlocks(cdsId, geneId, blockListextend, geneExon, cdsExon, intronList, geneLen,
                                         cdsLen, cdsSeq, geneSeq)
        
        listOfBlocExtendedExcluded = []
        if len(blockListextended) != 0:

            if len(blockListextended) >= 2:

                for i in range(0, len(blockListextended) - 1):
                    if len(blockListextended[i]) != 0 and len(blockListextended[i + 1]) != 0:
                        choiceNum = chooseBestIdentity(blockListextended[i], blockListextended[i + 1])
                        if choiceNum == 1:
                            listOfBlocExtendedExcluded.append(blockListextended[i + 1])
                        elif choiceNum == 2:
                            listOfBlocExtendedExcluded.append(blockListextended[i])
                        else:
                            pass
                    else:
                        pass

            listOfBlocExtendedNotExcluded = []
           
            for i in blockListextended:

                if len(i) != 0 and (i not in listOfBlocExtendedExcluded):
                    if len(i[0]) > 4:

                        listOfBlocExtendedNotExcluded.append(i[0][:-1])
                    else:
                        listOfBlocExtendedNotExcluded.append(i[0])
                else:
                    pass


            if len(listOfBlocExtendedNotExcluded) != 0:
                blockList = blockListLocal + listOfBlocExtendedNotExcluded
                blockList = order(blockList)
               
                blockList = concatenate_exons(blockList, cdsSeq, geneSeq)

            else:

                blockList = order(blockListLocal)
                blockList = concatenate_exons(blockList, cdsSeq, geneSeq)
        else:

            blockList = order(blockListLocal)
            blockList = concatenate_exons(blockList, cdsSeq, geneSeq)
    else:

        blockList = order(blockListLocal)
        blockList = concatenate_exons(blockList, cdsSeq, geneSeq)

   
    return blockList

def AlignOnlyRestOfExons(cdsId, geneId, blockListLocal, cdsExon, geneLen, cdsLen, geneExon, intronList, cdsSeq, geneSeq ):
    """
    This function do the PD alignment after founding the blocks of whole exon to align

    Parameters
    ----------

    geneId: string
        gene Id
    cdsId: string
        cds Id

    blockListLocal:
    cdsExon :
        
    geneLen: int
		gene length
    cdsLen: int
		cds length

    cdsSeq: string
	sequence of CDS 
    geneSeq: string
	sequence of gene
    geneExon:
    intronList:
   
    Returns
    -------
    blockList:list
        list of blocks
    """

    blockListDP = []
    blockListPDFound = []

   

    blockListDP = extarctBlockListToDP(cdsId, geneId, blockListLocal, cdsExon, geneLen, cdsLen)
    listToremove = []
    listToAdd = []

    if  len(blockListDP) !=0:

        if len(blockListDP) >= 2:
            for i in range (0, len(blockListDP)-1):

                bloc1= blockListDP[i]
                bloc2 = blockListDP[i+1]

                if bloc1[2] == bloc2[2] and bloc1[3] == bloc2[3]:


                    sameExon= True
                else:
                    sameExon = False
                if sameExon:
                    listToremove.append(i)
                    listToremove.append( i + 1)
                else:
                    listToremove.append('new')
            listOfElementToRemove = []
            if len(listToremove) != 0:

                listToAdd= filter(cdsId, geneId, blockListDP, listToremove)
                blockListDP= organize(blockListDP, listToremove, listToAdd)

            else:
                pass
            blockListDP = order(blockListDP)

        else:
            pass

        
        blockListPDFound = InterBlocksPD(blockListLocal, blockListDP, geneExon, cdsExon, intronList, cdsId, geneId, geneLen,
                                    cdsLen, cdsSeq, geneSeq)
    else:
       
        blockListPDFound = blockListLocal

    blockList = blockListPDFound
   
    return blockList


def filter( cdsId, geneId, blockListDP, listToremove):
    """
    This function

    Parameters
    ----------

    geneId: string
        gene Id
    cdsId: string
        cds Id
    blockListDP: list
    listToremove : list
    Returns
    -------
    listToAdd:list
        list of blocks
    """

    listToAdd = []

    if len(listToremove) != 0:

            i = 0
            while i < len(listToremove):

                if listToremove[i] == 'new':
                    while i < len(listToremove) and listToremove[i] == 'new':
                        
                        if i >= len(listToremove):
                            break
                    	i += 1

                if i >= len(listToremove):
                    break
               
                beginCds = blockListDP[listToremove[i]][0]
                beginGene = blockListDP[listToremove[i]][2]
                endGene = blockListDP[listToremove[i]][3]
               
                while i < len(listToremove) - 1 and listToremove[i] != 'new':
                    i += 1
              
                if i == len(listToremove) - 1 and listToremove[i] != 'new':

                    endCds = blockListDP[listToremove[i]][1]
                    listToAdd.append([beginCds, endCds, beginGene, endGene])
                    break
                elif listToremove[i] == 'new':
                    
                    endCds = blockListDP[listToremove[i - 1]][1]
                    listToAdd.append([beginCds, endCds, beginGene, endGene])
                    i += 1

    return listToAdd


def organize(blockListDP, listToremove, listToAdd):
    
    elementToRemove = []
    for i in listToremove:

        if type(i) == int and blockListDP[i] in blockListDP:
            elementToRemove.append(blockListDP[i])
        else:
            pass
    newelementToRemove = []
    for elt in elementToRemove:
        if elt not in newelementToRemove:
            newelementToRemove.append(elt)
    for bloc in newelementToRemove:
        blockListDP.remove(bloc)
    for bloc in listToAdd:
        blockListDP.append(bloc)
    

    return blockListDP



def ExtendAndAlignRestOfExons(cdsId, geneId, blockListLocal, cdsLen, cdsSeq, cdsExonForCds, cdsExon, geneLen, geneSeq, geneExon,
                              intronList):
    """
    This function extract blocks (rest of exons segments that weren't aligned by the local alignement ) and and do the extension 
			and the PD alignment

    Parameters
    ----------

    geneId: string
		gene Id
    cdsId: string
		cds Id
    blockListLocal: list
        list of blocks found by the local alignment
    cdsExon: list
		list of exons of cds
    geneLen: int
		gene length
    cdsLen: int
		cds length
    cdsExonForCds:
    cdsSeq:
    geneSeq:
    geneExon:
    intronList:


    Returns
    -------
    blockListToextend:list
        list of blocks extracted and that will be extended
    
    """

    blockListextended = []

    blockListPDFound = []


    blockListextended = OnlyExtend(cdsId, geneId, blockListLocal, cdsLen, cdsSeq, cdsExonForCds, geneLen, geneSeq, geneExon,
                           intronList)
   
    blockListPDFound = AlignOnlyRestOfExons(cdsId, geneId,  blockListextended , cdsExon, geneLen, cdsLen, geneExon, intronList,
                                     cdsSeq, geneSeq)

    return blockListPDFound


def extractBlockListToextend(cdsId, geneId,blockListLocal, cdsExon,geneLen,cdsLen):
    """
    This function extract blocks (rest of exons segments that weren't aligned by the local alignement ) and that will be extend

    Parameters
    ----------

    geneId: string
		gene Id
    cdsId: string
		cds Id
    blockListLocal: list
        list of blocks found by the local alignment
    cdsExon: list
		list of exons of cds
    geneLen: int
		gene length
    cdsLen: int
		cds length
    Returns
    -------
    blockListToextend:list
        list of blocks extracted and that will be extended
   
    """

    listOfBlock = blockListLocal
    blockListToextend = []
    for exon in cdsExon:

        for bloc in range(0, len(listOfBlock)):
           
            if exon[0] >= listOfBlock[bloc][1] or exon[1] <= listOfBlock[bloc][0]:
                pass
            elif exon[0] == listOfBlock[bloc][0]:
                if exon[1] == listOfBlock[bloc][1]:
                    pass
                else:

                    if ((bloc + 1) < len(listOfBlock)):
                        if listOfBlock[bloc + 1][0] >= exon[1]:
                            blockListToextend.append(
                                [listOfBlock[bloc][1], listOfBlock[bloc + 1][0], listOfBlock[bloc][3],
                                 listOfBlock[bloc + 1][2], 'endExon'])
                        else:

                            blockListToextend.append(
                                [listOfBlock[bloc][1], listOfBlock[bloc + 1][0], listOfBlock[bloc][3],
                                 listOfBlock[bloc + 1][2], 'middle'])
                    else:
                        blockListToextend.append([listOfBlock[bloc][1], exon[1], listOfBlock[bloc][3], geneLen, 'endExon'])
            else:

                if exon[1] == listOfBlock[bloc][1]:

                    if bloc != 0:
                        if listOfBlock[bloc - 1][0] >= exon[0]:
                            blockListToextend.append(
                                [listOfBlock[bloc - 1][1], listOfBlock[bloc][0], listOfBlock[bloc - 1][3],
                                 listOfBlock[bloc][2], 'middle'])
                        else:
                            blockListToextend.append(
                                [listOfBlock[bloc - 1][1], listOfBlock[bloc][0], listOfBlock[bloc - 1][3],
                                 listOfBlock[bloc][2], 'beginExon'])
                    else:
                        blockListToextend.append([exon[0], listOfBlock[bloc][0], 0, listOfBlock[bloc][2], 'beginExon'])
                else:

                    if bloc != 0:

                        if listOfBlock[bloc - 1][0] >= exon[0]:
                            blockListToextend.append([listOfBlock[bloc - 1][1], listOfBlock[bloc][0], listOfBlock[bloc - 1][3],
                                                      listOfBlock[bloc][2], 'middle'])
                        else:
                            blockListToextend.append([listOfBlock[bloc - 1][1], listOfBlock[bloc][0], listOfBlock[bloc - 1][3],
                                                      listOfBlock[bloc][2], 'beginExon'])

                    else:
                        blockListToextend.append([exon[0], listOfBlock[bloc][0], 0, listOfBlock[bloc][2], 'beginExon'])

                    if ((bloc + 1) < len(listOfBlock)):
                        if listOfBlock[bloc + 1][0] >= exon[1]:
                            blockListToextend.append([listOfBlock[bloc][1], listOfBlock[bloc + 1][0], listOfBlock[bloc][3],
                                                      listOfBlock[bloc + 1][2], 'endExon'])
                        else:

                            blockListToextend.append([listOfBlock[bloc][1], listOfBlock[bloc + 1][0], listOfBlock[bloc][3],
                                                  listOfBlock[bloc + 1][2], 'middle'])
                    else:
                        blockListToextend.append([listOfBlock[bloc][1], exon[1], listOfBlock[bloc][3], geneLen, 'endExon'])

    newblockListToextend = []
    for i in blockListToextend:
        if i not in newblockListToextend:
            newblockListToextend.append(i)
   
    return newblockListToextend

def chooseBestIdentity(blocks1, blocks2):
    """
    This function chooses in case of overlapping of two block, the bloc that has higher identity ,

    Parameters
    ----------
    blocks1: list
        bloc1 with its identity
    blocks2: list
       bloc2 with its identity

    Returns
    -------

    int: 1 if bloc1, 2 if bloc2 an3 if both and no overlapping
    """

    [bloc1, identity1] = blocks1
    [bloc2, identity2] = blocks2

    if bloc1[3]>= bloc2[2] :
        if identity1 >= identity2:
            return   1
        else:
            return  2
    else:
        return 3


def extendBlocks(cdsId, geneId,blockListToextend, geneExon, cdsExon, intronList, geneLen, cdsLen,
                                        cdsSeq, geneSeq):
    """
    This function try to extend the segment of exon

    Parameters
    ----------

   geneId: string
		gene id
   cdsId: string
	cds id

   geneSeq: string
	gene sequence
   cdsSeq: string
	cds sequence
   cdsExon: list
	list of exons of cds
    blockListToextend: list
        bloc that will be extended
    intronList:dictionary
               dictionary of gene intron, the keys are gene Id and the values are the blocks interval of intron gene
    geneExon: list
        list of exon of gene

    geneLen: int
		gene length
	cdsLen: int
		cds length
    Returns
    -------
    blockListextend: list
              the extended bloc
    """




    blockListextend=[]
    for bloc in blockListToextend:


        [cdsBegin, cdsEnd, genBegin, geneEnd, position] = bloc
        if position ==  'endExon':
            blockListextend.append(extendLeftToRight( bloc, cdsId, geneId,  cdsSeq, geneSeq))
        elif position ==  'beginExon':
            blockListextend.append(extendRightToLeft( bloc, cdsId, geneId,  cdsSeq, geneSeq))
	else:
		blockListextend.append(extendLeftToRight_Middle( bloc, cdsId, geneId,  cdsSeq, geneSeq))
		
    return blockListextend
def extendLeftToRight_Middle(bloc, cdsId, geneId,cdsSeq, geneSeq):
    """
    This function try to extend the middle segment of exon doing steps from left to rigth

    Parameters
    ----------

    bloc: list
        bloc to extend
   geneId: string
        gene id
    cdsId: string
        cds id

    geneSeq: string
        gene sequence
    cdsSeq: string
        cds sequence
    
    Returns
    -------
    blockListextend: list
              the extended bloc
    """
    blockListextend = []
    [cdsBegin, cdsEnd, genBegin, geneEnd, position] = bloc
    remain = abs(geneEnd - genBegin) / abs(cdsEnd - cdsBegin)
    if 0.5 <= remain and remain <= 2:

       blockListextend=[bloc, CUTOFF_IDENTITY]

    elif remain > 2:

        differenceSequence = abs(geneEnd - genBegin) - abs(cdsEnd - cdsBegin)

        
        if differenceSequence >= ALPHA:
            blockListextend = computeBestMatch(cdsId, geneId, ALPHA, cdsSeq, cdsBegin, cdsEnd, geneSeq, genBegin,
                                               geneEnd)
        else:
            blockListextend = computeBestMatch(cdsId, geneId, differenceSequence, cdsSeq, cdsBegin, cdsEnd, geneSeq,
                                               genBegin, geneEnd)
    else:
	pass
    return blockListextend

def extendLeftToRight(bloc, cdsId, geneId,cdsSeq, geneSeq):
    """
    This function try to extend the segment of exon doing steps from left to rigth

    Parameters
    ----------

    bloc: list
        bloc to extend
   geneId: string
        gene id
    cdsId: string
        cds id

    geneSeq: string
        gene sequence
    cdsSeq: string
        cds sequence
    
    Returns
    -------
    blockListextend: list
              the extended bloc
    """
    blockListextend = []
    [cdsBegin, cdsEnd, genBegin, geneEnd, position] = bloc
    if abs(geneEnd - genBegin) < abs(cdsEnd - cdsBegin):

       pass
    
    elif abs(geneEnd - genBegin) == abs(cdsEnd - cdsBegin):
        cdsseqtemp = cdsSeq[cdsBegin:cdsEnd]
        geneseqtemp = geneSeq[genBegin:geneEnd]
       
        identityAlignment = computeAlignmentPercentIdentity(cdsseqtemp, geneseqtemp)

        if identityAlignment >= CUTOFF_IDENTITY:
            blockListextend=[bloc, identityAlignment]

        else:

            pass
    
    else:

        differenceSequence = abs(geneEnd - genBegin) - abs(cdsEnd - cdsBegin)

        
        if differenceSequence >= ALPHA:
            blockListextend = computeBestMatch(cdsId, geneId, ALPHA, cdsSeq, cdsBegin, cdsEnd, geneSeq, genBegin,
                                               geneEnd)
        else:
            blockListextend = computeBestMatch(cdsId, geneId, differenceSequence, cdsSeq, cdsBegin, cdsEnd, geneSeq,
                                               genBegin, geneEnd)
    return blockListextend

def extendRightToLeft( bloc, cdsId, geneId, cdsSeq, geneSeq):
    """
     This function try to extend the segment of exon doing steps in the right side of gene

     Parameters
     ----------

     bloc: list
         bloc to extend
    geneId: string
         gene id
     cdsId: string
         cds id

     geneSeq: string
         gene sequence
     cdsSeq: string
         cds sequence

     Returns
     -------
     blockListextend: list
               the extended bloc
     """
    blockListextend = []
    [cdsBegin, cdsEnd, genBegin, geneEnd, position] = bloc
   
    if abs(geneEnd - genBegin) < abs(cdsEnd - cdsBegin):

        pass
    elif abs(geneEnd - genBegin) == abs(cdsEnd - cdsBegin):
        cdsseqtemp = cdsSeq[cdsBegin:cdsEnd]
        geneseqtemp = geneSeq[genBegin:geneEnd]
        
        identityAlignment = computeAlignmentPercentIdentity(cdsseqtemp, geneseqtemp)

        if identityAlignment >= CUTOFF_IDENTITY:
            blockListextend=[bloc, identityAlignment]

        else:

            pass
    else:


        differenceSequence = abs(geneEnd - genBegin) - abs(cdsEnd - cdsBegin)
       
        if differenceSequence >= ALPHA:
            
            genBegin = geneEnd- (len(cdsSeq[cdsBegin:cdsEnd]) * ALPHA) +ALPHA

            blockListextend = computeBestMatch(cdsId, geneId, ALPHA, cdsSeq, cdsBegin, cdsEnd, geneSeq, genBegin, geneEnd)
        else:
            genBegin = geneEnd - (len(cdsSeq[cdsBegin:cdsEnd]) * differenceSequence)+differenceSequence

            blockListextend = computeBestMatch(cdsId, geneId, differenceSequence, cdsSeq, cdsBegin, cdsEnd, geneSeq,
                                               genBegin, geneEnd)

    return blockListextend

def computeBestMatch(cdsId, geneId, ALPHAInput, cdsSeq, cdsBegin, cdsEnd, geneSeq, genBegin, geneEnd):
    """
     This function compute the best segment and  try to extend the segment taking account the best step

     Parameters
     ----------


    geneId: string
         gene id
     cdsId: string
         cds id

     geneSeq: string
         gene sequence
     cdsSeq: string
         cds sequence
    cdsBegin: int
        cds begin
    cdsEnd:int
        cds end
    genBegin:int
        gene begin
    geneEnd: int
        gene end

    ALPHAInput: int
        max of step
     Returns
     -------
     blockListextend: list
               the extended bloc
     """
    blockListextend = []
    identityAlignmentCut = 0

    
    if ALPHAInput%3==0:
        ALPHA_ = ALPHAInput

    elif (ALPHAInput-1)%3 == 0:
                ALPHA_ = ALPHAInput -1
    elif (ALPHAInput-2)%3 == 0:
                ALPHA_ = ALPHAInput -2


    blockListextend= chooseBestBlock(3, cdsBegin, cdsEnd, genBegin, geneEnd , cdsSeq, geneSeq, ALPHA_)
    return blockListextend


def chooseBestBlock(step, cdsBegin, cdsEnd, genBegin, geneEnd , cdsSeq, geneSeq, ALPHA_):
	identityAlignmentCut=0
	blockListextend=[]
        k=0
	for i in range(cdsBegin, cdsEnd, step):
		cdsseqtemp = cdsSeq[i:cdsEnd]
		for j in range (genBegin,geneEnd-len(cdsseqtemp), 3 ):
			k += 1
			if k > ALPHA_:
				break
			
			geneseqtemp = geneSeq[j:j+ len(cdsseqtemp)]
	       		
		
			identityAlignment = computeAlignmentPercentIdentity(cdsseqtemp, geneseqtemp)
		       
			if identityAlignment > identityAlignmentCut:
			    identityAlignmentCut = identityAlignment
			    bestStep = i
		   
		if identityAlignmentCut >= CUTOFF_IDENTITY:
			blockListextend=[[i, cdsEnd, j,j+ len(cdsseqtemp)], identityAlignmentCut]
			
			break
	return blockListextend

def computeAlignmentPercentIdentity(sequence1, sequence2):
    """
    This function aim to return the percentage of identity of 2 sequences
    sequence1: string
        sequence 1
    sequence2: string
        sequence 2

     Returns
     -------
     aln_identity: float
        % of identity between the 2 sequences
    """

    aln_identity = 0.0
    match = 0
    length = len(sequence1)

    for i in range(length):

        if (sequence1[i] == sequence2[i]):
            match += 1
    aln_identity = 100.0 * match / length
    return aln_identity


def extarctBlockListToDP(cdsId, geneId,blockListLocal, dictcdsExon,geneLen,cdsLen):
    """
    This function extract blocks (rest of exons that weren't aligned by the local alignement ) and that will be forced to be aligned

    Parameters
    ----------
    geneId: string
         gene id
     cdsId: string
         cds id

    blockListLocal: list
        list of blocks found by the local alignment
    cdsExon: list
		list of exons of cds
    geneLen: int
		gene length
    cdsLen: int
		cds length
    Returns
    -------
    blockListPD:list
        list of blocks extracted and that will be extended


    """
    blockListPD = []
    listExonLosses = []
    cdsExon =     dictcdsExon [cdsId]

    if len(blockListLocal) == 0:
        for exon in cdsExon:
            blockListPD.append([exon[0], exon[1], 0, geneLen])
        return blockListPD
    else:
        for exon in cdsExon:
            for bloc in range(0, len(blockListLocal)):

                beginBlocCds = blockListLocal[bloc][0]
                endBlocCds = blockListLocal[bloc][1]
                if exon[0] <= beginBlocCds and endBlocCds <= exon[1]:
                    break
                if bloc == len(blockListLocal) - 1:
                    listExonLosses.append(exon)

        if len(listExonLosses) == 0:
            return []
        else:

            for exon in listExonLosses:

                precedent = [0, 0, 0, 0]
                next = [exon[0], exon[1], geneLen, geneLen]
                for bloc in range(0, len(blockListLocal)):

                    if exon[0] >= blockListLocal[bloc][1]:
                        precedent = blockListLocal[bloc]

                    if exon[1] <= blockListLocal[bloc][0]:
                        next = blockListLocal[bloc]
                        break

                blockListPD.append([exon[0], exon[1], precedent[3], next[2]])

        return blockListPD


def analyse_result(blockList, cdsId, cdsSequence, geneId, geneSequence, sourceData, cds2GeneExon, cdsExon):
    """
    This function analyses results of blocks found

    Parameters
    ----------


    blockList: list
    list of blocks found by the method of alignment
    geneId: string
            gene Id
    geneSequence: string
            gene sequence
    cds: string
        CDS Id
    cdsSequence: string
            CDS sequence

    cdsExon: dictionary
       dictionay with CDS id as key and exons interval list as value
    sourceData: list
       list that contain information about CDS like id CDS, CDS sequence and id of its gene

    cds2GeneExon: dictionary
       dictionay with cds id as key and exons blocks list as values


    Returns
    -------
    status: int

    splicing_sites: list
    list ofsplicing sites
    targetcds: string

    texisting: int

    """

    status = -1
    texisting = -1

   
    # create the target of the cds
    targetcds = create_target_cds(blockList, geneSequence)

    # compute the splicing sites
    splicing_sites = compute_splicing_sites(blockList, geneSequence)

    # create the target of exons
    targetexon = create_target_exon(blockList)

    
    if (len(blockList) > 0 and is_continuous_onCDS(blockList)):
        ccdsexon = cdsExon[cdsId]
        existing = False
        i = 0
        while(i < len(sourceData) and not existing):
            tcds = sourceData[i]
            tcdsid,tcdsseq,tgeneid,null = tcds
            ccdsexon = cdsExon[cdsId]
            tcdsexon = cdsExon[tcdsid]
            geneexon = cds2GeneExon[tcdsid]
            acceptor = True
            donor = True
            frame = True
            if(tgeneid == geneId and len(ccdsexon) == len(tcdsexon) == len(targetexon) == len(geneexon)):
                if(targetexon[0][1] != geneexon[0][1]):
                    donor = False
                for j in range(1,len(targetexon)-1):
                    if(targetexon[j][0] != geneexon[j][0]):
                        acceptor = False
                    if(targetexon[j][1] != geneexon[j][1]):
                        donor = False
                if(targetexon[-1][0] != geneexon[-1][0]):
                    acceptor = False
                for j in range(len(ccdsexon)):
                    sizeccds = ccdsexon[j][1]- ccdsexon[j][0]
                    sizetcds = tcdsexon[j][1]- tcdsexon[j][0]
                    if(sizeccds - sizetcds) % 3 != 0:
                        frame = False
            else:
                frame = False

            if(acceptor == True and donor == True and frame == True):
                existing = True
                status = STATUS_EXISTING_PROTEIN
                texisting = i
            i+=1
        if(not existing):
            if len(targetexon) == len(cdsExon[cdsId]):
                if(is_protein(targetcds)):
                    status = STATUS_PREDICTED_PROTEIN
                else:
                    status = STATUS_PREDICTED_CDS
            else:
                status = STATUS_COMPLETE_CDS_DIFFERENT_STRUCTURE
    else:
        status = STATUS_PARTIAL_CDS


    return status, splicing_sites,targetcds,texisting




##############################
### ORTHOLOGY COMPUTATION ####
##############################

def compute_orthology_matrix(sourcedata,targetdata,comparisonresults):
    """


    Parameters
    ----------
    sourceData: list
    	list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
    	list that contain list of all informations about all genes:  id and sequence of each gene

    comparisonresults: list
	blocs result from used alignment method 


    Returns
    -------
    orthologymatrix:matrix


    """
    orthologymatrix = []
    for cds in sourcedata:
        orthologymatrix.append([])
        for cds in sourcedata:
            orthologymatrix[-1].append(0)

    i = 0
    for gene in targetdata:
        geneid,geneseq = gene
        j = 0
        for cds in sourcedata:
            cdsid,cdsseq,cdsgeneid, null = cds
            status, blocklist,splicing_sites,targetcds, texisting = comparisonresults[i][j]
            if(status == STATUS_EXISTING_PROTEIN):
                orthologymatrix[j][texisting] = 1
                orthologymatrix[texisting][j] = 1
            j += 1
        i += 1
    return orthologymatrix

def computeOrthology(sourcedata,targetdata,comparisonresults):
    """


     Parameters
     ----------
     sourceData: list
    	list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
    	list that contain list of all informations about all genes:  id and sequence of each gene

    comparisonresults: list
	blocs result from used alignment method 




     Returns
     -------
     grouplist:list 
     """
    print "Computing orthology groups"
    grouplist = []
    
    orthologymatrix = compute_orthology_matrix(sourcedata,targetdata,comparisonresults)
    for gene1 in range(len(orthologymatrix)):
        assigned_to_group = False
        for i in range(len(grouplist)):
            group = grouplist[i]
            group_assigned = False
            for gene2 in group:
                if(orthologymatrix[gene1][gene2] == 1):
                    group_assigned = True
            if(group_assigned):
                assigned_to_group = True
                grouplist[i].append(gene1)
    
        if(not assigned_to_group):
            grouplist.append([gene1])

    grouplist = merge_orthology_groups(grouplist)
    return grouplist

def merge_orthology_groups(grouplist):
    """


     Parameters
     ----------
     grouplist:

     Returns
     -------
     grouplist:
     """
    groups_to_delete = []
    for i in range(len(grouplist)):
        for j in range(i+1,len(grouplist)):
            if(len(list(set(grouplist[i]) & set(grouplist[j]))) != 0):
                grouplist[i] += grouplist[j]
                groups_to_delete.append(grouplist[j])

    newgrouplist=[]
    for sublist in groups_to_delete:
        if sublist not in newgrouplist:
            newgrouplist.append(sublist)
    
    for group in newgrouplist:
        grouplist.remove(group)
    
    return grouplist


def completeOrthology(sourcedata,targetdata,comparisonresults,orthologygroups):
    """


     Parameters
     ----------

    sourceData: list
    	list that contain information about CDS like id CDS, CDS sequence and id of its gene
    targetData: list
    	list that contain list of all informations about all genes:  id and sequence of each gene

    comparisonresults: list
	blocs result from used alignment method 

     orthologygroups:


     Returns
     -------
     sourcedata:
     orthologygroups:

     """
    print "Adding predicted protein-coding sequence"
    cds2id = []
    cds2geneid = []
    geneidlist = []
    
    nbcdsinitial = len(sourcedata)
    for cds in sourcedata:
        cdsid,cdsseq,cdsgeneid, null = cds
        cds2id.append(cdsid)
        cds2geneid.append(cdsgeneid)

    for k in range(len(orthologygroups)):
        group2geneid = []
        for cds in orthologygroups[k]:
            group2geneid.append(cds2geneid[cds])
        
        for i in range(len(targetdata)):
            gene = targetdata[i]
            geneid,geneseq = gene
            if(geneid not in group2geneid):
                #predicted = False
                for cds in orthologygroups[k]:
                    if(cds < nbcdsinitial):
                        cdsid,null,null,null = sourcedata[cds]
                        status, blocklist,splicing_sites,targetcds, texisting = comparisonresults[i][cds]
                        all_valid_sites = all_valid_splicing_sites(splicing_sites)
                        if((status == STATUS_PREDICTED_PROTEIN or status == STATUS_COMPLETE_CDS) and all_valid_sites):
                            prediction = ""
                            if(status == STATUS_PREDICTED_PROTEIN):
                                prediction = "Protein"
                            else:
                                prediction = "CDS"
                                
                            pid =  "Predicted_"+prediction+"_Group"+str(k)+"_"+cdsid+"_to_"+geneid
                            pseq = targetcds
                            pgeneid = geneid
                            new_prediction = True
                            j = len(sourcedata)-1
                            while(j >= 0 and sourcedata[j][2] == pgeneid):
                                new_prediction = new_prediction and (pseq != sourcedata[j][1])
                                j -= 1
                            if(new_prediction):
                                sourcedata.append([pid,pseq,pgeneid,[]])#TO COMPUTE
                                orthologygroups[k].append(len(sourcedata)-1)
                                
    return sourcedata,orthologygroups


