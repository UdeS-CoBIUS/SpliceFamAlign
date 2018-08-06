#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``alignment.py`` **module description**:


This module defines functions required at each step of the comparison.

.. moduleauthor:: Aïda Ouangraoua, Safa Jammali and Jean-David Aguilar
Université de Sherbrooke Canada
Laboratoty CoBiUS

2017-2018

"""

import os
import csv
import argparse
from Bio import pairwise2
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

from blocklistManagement import *

ALPHA_TRIM=5
MIN_IDENTITY_ALN = 0.5

MIN_EVALUE = 1.0/10000000 # 10-7
MAX_EVALUE = 0.1 # 10-1
DELTA_EVALUE = 100
MIN_CDSLEN_FOR_ALN = 15
MINLEN_FOR_JOINING = 6
MINIMAL_LENGTH_ALN_EXACT=25
MIN_BLOCK_EXTENSION_LEN = 10
MIN_BLOCK_EXTENSION_IDTY_DECREASE = 0.2

MAX_SHIFT = 30 # _FOR_SPLICE_JUNCTION_CORRECTION
MAX_GAP = 3 # _FOR_SPLICE_JUNCTION_CORRECTION

MAX_CDS_SEPARATION_FOR_CONCAT = 100
MAX_SEPARATION_DIFF_FOR_CONCAT = 60
MIN_VALIDBLOCK_LEN = 50
BORDER_LENGTH = 80
LENGTH_ALN_LINE = 80
GLOBAL_LIMIT = 10**6

def launch_splign(cdsseqfilename, geneseqfilename, cdsid, geneid):
    """
    This function launchs splign

    Parameters
    ----------
    cdsseqfilename:
    geneseqfilename:
    cdsid:
    geneid:

    Returns
    -------
    splignoutputname:
    """

    splignoutputname = os.getcwd()+'/src/results/splign_results/CDS_' + cdsid+'VS_Gene_' + geneid+".splign"

    os.system(os.getcwd()+"/src/splign -query "+ cdsseqfilename + " -subj "+geneseqfilename+" -disc -aln "+splignoutputname+ " -min_compartment_idty 0.00 -min_exon_idty 0.00")

    return splignoutputname

def launch_splign_default(cdsseqfilename, geneseqfilename, cdsid, geneid):
    """
    This function launch splign with default parameters

    Parameters
    ----------
    cdsseqfilename:
    geneseqfilename:
    cdsid:
    geneid:

    Returns
    -------
    splignoutputname:
    """
    splignoutputname = os.getcwd()+'/src/results/splign_results/CDS_' + cdsid+'VS_Gene_' + geneid+".splign"
    
    os.system(os.getcwd()+"/src/splign -query "+ cdsseqfilename + " -subj "+geneseqfilename+" -disc -aln "+splignoutputname)

    return splignoutputname

def parse_splign_output(splignoutputname):
    """
    This function launch splign

    Parameters
    ----------
    splignoutputname:

    Returns
    -------
    blocklist
    """
    file=open(splignoutputname, "r")
    blocklist = []
    for line in file.readlines(): 
        if 'No alignment was found' in line or ">-1" in line:
            break
        else:
            if("Exon" in line):
                location = line.split("(")[1].split(")")[0]
                identity= float(line.split(" ")[-1])
               
                querylocation,subjectlocation = location.split(",")
                querylocation = [int(x) for x in querylocation.split("-")]
                querylocation[0] -= 1
                subjectlocation = [int(x) for x in subjectlocation.split("-")]
                subjectlocation[0] -= 1
                block = querylocation + subjectlocation
                if identity >= 1:
                    blocklist.append(block)
    file.close()
    os.remove(splignoutputname)
    return blocklist

def BlastLocal(cdsSeqFile, cdsSeq, geneSeqFile,geneSeq,  cdsId, cdsGeneId, evalue ):
    """
    This function launchs blast in order to construct the structure of the input 

    Parameters
    ----------
    cdsSeqFile:
    cdsSeq:
    geneSeqFile:
    geneSeq:
    cdsId:
    cdsGeneId:
    evalue:

    Returns
    -------
    blockList_:
    blockList
    """
    listHits=[]
    blockList=[]
    blockList_=[]
    blocklistout= launch_tblastx(cdsSeqFile, geneSeqFile,evalue, [0,0,0,0],[0,0,0,0],cdsId, cdsGeneId)
    listHits= readBlastOut(blocklistout, 0,0)
    listHitsExonInOrder =sorted(listHits, key=lambda colonnes: colonnes[4])
    blockList= constructExonStructure(cdsSeq,geneSeq,listHitsExonInOrder )
    blockList=order( blockList)
    for i in range (0, len(blockList)):
	blockList_.append(blockList[i][0:4])
   
    os.remove(blocklistout)
    return blockList_, blockList

def readBlastOut(blastOutput, blockQueryEnd, blockSubjectEnd):
	"""
	This function reads blast output file and extract the hits 

	Parameters
	----------
	bastOutput:
	blockQueryEnd:
	blockSubjectEnd

	Returns
	-------
	listHits:
	"""
	listHits=[]
	with open(blastOutput, 'rb') as csvfile:
		spamreader = csv.reader(csvfile,delimiter = ' ', quotechar = ' ')

		for row in spamreader:


			if row[0][0] is not '#':
				# construct a hit

				identity= float(row[0].split('\t')[2])
				beginNewHitCDS= int(row[0].split('\t')[6])+int(blockQueryEnd) - 1
				endNewHitCDS = int(row[0].split('\t')[7]) + int(blockQueryEnd)
				beginNewHitgene = int(row[0].split('\t')[8])+ int(blockSubjectEnd) - 1
				endNewHitgene = int(row[0].split('\t')[9]) + int(blockSubjectEnd)
				evalueNewHit = float(row[0].split('\t')[10])
				if (beginNewHitCDS < endNewHitCDS and beginNewHitgene < endNewHitgene and evalueNewHit <= MIN_EVALUE):
					listHits.append([beginNewHitCDS,endNewHitCDS, beginNewHitgene, endNewHitgene, evalueNewHit, identity])
	return listHits

 

   
def constructExonStructure( cdsSeq,geneSeq, listHits ):
        """
	This function treats the hits in order to construct the structure of inputs 

	Parameters
	----------
	cdsSeq:
        geneSeq:
        listHits

	Returns
	-------
	listHits:
	"""
	blocks=[]
	countHit = 0
	while not cover_wholeCDS(blocks,len(cdsSeq)):
		
		i=0
		while i < len(listHits):
			
			listRemove=[]
			
			currentHit= treatHit( listHits[i], cdsSeq,geneSeq)
                        
			if len(currentHit)!=0:
				if len(blocks)>0:
					
					for j in range (0,len( blocks)):
						
						if blocks[j][0] <= currentHit[0] and currentHit[0] <= blocks[j][1] \
						and blocks[j][0] <= currentHit[1] and currentHit[1] <= blocks[j][1] \
						and blocks[j][2] <= currentHit[2] and currentHit[2] <= blocks[j][3] \
						and blocks[j][2] <= currentHit[3] and currentHit[3] <= blocks[j][3]:
							
							break
						# new hit and previous hit  are different
						elif (blocks[j][0] >= currentHit[1] and blocks[j][2] > currentHit[3]) \
								or (currentHit[0] >= blocks[j][1] and currentHit[2] > blocks[j][3]):
							
							if j == len(blocks)-1:
								
								blocks.append(currentHit)
													
							
							else:
								pass
			


						# new Hit and previous are overlaping, right block
						elif (min(currentHit[0], blocks[j][0]) == blocks[j][0] and min(currentHit[2], blocks[j][2]) == blocks[j][2]) :							
							if  currentHit[0]< blocks[j][1] and blocks[j][1]< currentHit[1] and currentHit[2]< blocks[j][3] and blocks[j][3]< currentHit[3]:


								
								currenthit= [blocks[j][0],currentHit[1] , blocks[j][2], currentHit[3],min(currentHit[4], blocks[j][4]), 100.00]
								listRemove.append(blocks[j])
								
								if j == len(blocks)-1:
									blocks.append(currentHit)
									break
							
							#trim based on which respect AG GT else keep best eval
							elif currentHit[0]< blocks[j][1] and blocks[j][1]< currentHit[1] and  blocks[j][3]< currentHit[2]: 
									#case=1
									case=trimAGGT(blocks[j],currentHit, cdsSeq,geneSeq )
									
									if case==1:
										currentHit = [currentHit[0]+(blocks[j][1]-currentHit[0]), currentHit[1], currentHit[2]+(blocks[j][1]-currentHit[0]), currentHit[3], currentHit[4], currentHit[5]]
										if j == len(blocks)-1:
											blocks.append(currentHit)
											break
								
									else:
										newBloc=[blocks[j][0], blocks[j][1]-(blocks[j][1]-currentHit[0]), blocks[j][2], blocks[j][3]-(blocks[j][1]-currentHit[0]), blocks[j][4], blocks[j][5]]
										blocks[blocks.index(blocks[j])] = newBloc
										
										if j == len(blocks)-1:
											blocks.append(currentHit)
											break
									

							#keep best eval
							else:
								
								if min(currentHit[4], blocks[j][4])== blocks[j][4]:
									break
								else:
									
									if j == len(blocks)-1:
										listRemove.append(blocks[j])
										blocks.append(currentHit)
										break
									else:
								
										listRemove.append(blocks[j])
								
						# new Hit and previous are overlaping, left block
						elif (min(currentHit[0], blocks[j][0])  == currentHit[0] and min(currentHit[2], blocks[j][2])  == currentHit[2]):
							
							#concat both
							if blocks[j][0]< currentHit[1] and currentHit[1]< blocks[j][1] and blocks[j][2]< currentHit[3] and currentHit[3]< blocks[j][3]:

								
								currentHit= [currentHit[0] , blocks[j][1], currentHit[2],  blocks[j][3], min(currentHit[4], blocks[j][4]), 100.00]
								
								listRemove.append(blocks[j])
								
								if j == len(blocks)-1:
									blocks.append(currentHit)
									break
							#trim based on which respect AG GT else keep best eval
							elif blocks[j][0]< currentHit[1] and currentHit[1]< blocks[j][1] and  currentHit[3]< blocks[j][2]: 

								case=trimAGGT(currentHit,blocks[j], cdsSeq,geneSeq )
								
								if case==1:
									newBloc= [blocks[j][0]+(currentHit[1]-blocks[j][0]), blocks[j][1], blocks[j][2]+(currentHit[1]-blocks[j][0]), blocks[j][3], blocks[j][4], blocks[j][5]]
									blocks[blocks.index(blocks[j])] = newBloc
									
									if j == len(blocks)-1:
										blocks.append(currentHit)
										break
									
								else:
									currentHit = [currentHit[0], currentHit[1]-(currentHit[1]-blocks[j][0]), currentHit[2], currentHit[3]-(currentHit[1]-blocks[j][0]), currentHit[4], currentHit[5]]
									if j == len(blocks)-1:
										blocks.append(currentHit)
										break
								
							#keep best eval
							else:
								
								if min(currentHit[4], blocks[j][4])== blocks[j][4]:
									break
								else:
									
								
									if j == len(blocks)-1:
										listRemove.append(blocks[j])
										blocks.append(currentHit)
										break
									else:
							
										listRemove.append(blocks[j])
								
			
						#keep best e val
						else:
							
							if min(currentHit[4], blocks[j][4])== blocks[j][4]:
									break
							else:
									
								
									if j == len(blocks)-1:
										listRemove.append(blocks[j])
										blocks.append(currentHit)
										break
									else:
							
										listRemove.append(blocks[j])
								

					if len(listRemove)>0:
						for element in listRemove:
							blocks.remove(element)	

				# add first hit after doing treatments: extend or trim based on percentIdentity
				else:
			
				
					blocks.append(currentHit)


				
			else:
				pass
				
			i= i+1
			
		
		
		if i >= len(listHits):
			break
		
	blocks = order(blocks)	
  	
	# if finish hits and not yet cover 100% then add remainnig exons
	if not cover_wholeCDS(blocks,len(cdsSeq)):
		blocks, remainingBlocs= extractandExtend(blocks, cdsSeq,geneSeq)
		

		 #verif not cds not covered di local PD alignment
		if not cover_wholeCDS(blocks,len(cdsSeq)):
			
			for j in remainingBlocs:
				blockToadd=localAlignDP(j, cdsSeq,geneSeq)
				for bl in blockToadd:
					blocks.append(bl)
			blocks=order(blocks)
			if not cover_wholeCDS(blocks,len(cdsSeq)):
				print 'Not Covered all CDS'
				exit(-1)
			else:
				print 'Covered all CDS'
		else:
				print 'Covered all CDS'
	else:
				print 'Covered all CDS'
		
	return blocks    
def localAlignDP(bloc, cdsSeq,geneSeq):
        """
	This function computes local alignment in order to complete the structure of inputs 

	Parameters
	----------
	cdsSeq:
        geneSeq:
        bloc

	Returns
	-------
	Foundblocs :
	"""
	se1=geneSeq[bloc[2]:bloc[31]]
	se2=cdsSeq[bloc[0]:bloc[1]]
	alignment =pairwise2.align.localxx(se1, se2)
	state= False
	begin=[]
	end=[]
	dicAlignmentPosandCdsGeneNt=  {}
	ntC =0
	ntG =0
	for i, j in enumerate (alignment[0][1]):

		if j!='-':
			ntC +=1
		if alignment[0][0][i]!='-':
			ntG +=1
		if  j==alignment[0][0][i]:
			
			if   not state:
				begin.append(i)
				dicAlignmentPosandCdsGeneNt[i]=(ntC-1, ntG-1)
				state= True
				if i==len(alignment[0][1])-1:
					end.append(i+1)
					dicAlignmentPosandCdsGeneNt[i+1]=(ntC, ntG)
			else:
							
				if i==len(alignment[0][1])-1:
					end.append(i+1)
					dicAlignmentPosandCdsGeneNt[i+1]=(ntC, ntG)
		else:
			
			if   not state:
			
				pass
			
			else:
				end.append(i)
				dicAlignmentPosandCdsGeneNt[i]=(ntC, ntG-1)
				state= False
	
	Foundblocs=[]
	if len(begin)>0:
		for i in range (0, len(begin)):
			Foundblocs.append([ dicAlignmentPosandCdsGeneNt[begin[i]][0]+bloc[0], dicAlignmentPosandCdsGeneNt[end[i]][0]+bloc[0], dicAlignmentPosandCdsGeneNt[begin[i]][1]+bloc[2], dicAlignmentPosandCdsGeneNt[end[i]][1] +bloc[2]])
		
	return Foundblocs 
	
	
def  treatHit(currentHit, cdsSeq,geneSeq):
	"""
	This function 
	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:

	Returns
	-------
	hit :
	"""
	percentIdentity=currentHit[5]
	if percentIdentity == 100.00:
			
			hit = extendHit(currentHit, cdsSeq,geneSeq)
			

	else: 
			
			hit= trimHit(ALPHA_TRIM, currentHit, cdsSeq,geneSeq)
			if len (hit)> 0:
				hit= extendHit(hit, cdsSeq,geneSeq)
			

	return hit

def trimHit(ALPHA_TRIM, currentHit, cdsSeq,geneSeq):
	"""
	This function 

	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:
	ALPHA_TRIM:
	Returns
	-------
	hit :
	"""
	percentIdentity= currentHit[5]
	
	if len(cdsSeq[currentHit[0]: currentHit[1]-ALPHA_TRIM])<= ALPHA_TRIM:
		
			return []
	else:
			
			hit= trimRigth(ALPHA_TRIM, currentHit, cdsSeq,geneSeq)
			if len(hit)==0:
				hit= trimLeft(ALPHA_TRIM, currentHit, cdsSeq,geneSeq)
				if len(hit)==0:
					hit= trimRightLeft(ALPHA_TRIM, currentHit, cdsSeq,geneSeq)
					if len(hit)==0:
						return []
					else: 
						return hit
				else: 
					return hit
			else: 
				return hit
			
			
			
def trimRigth(ALPHA_TRIM, currentHit, cdsSeq,geneSeq) :
	"""
	This function local alignment in order to complete the structure of inputs 

	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:
	ALPHA_TRIM:
	Returns
	-------
	hit :
	"""
	for i in range (currentHit[1]-ALPHA_TRIM, currentHit[0]+ALPHA_TRIM, -ALPHA_TRIM ):
		
		percentIdentity=compute_alignment_identity(cdsSeq[currentHit[0]: i],geneSeq[currentHit[2]: currentHit[3]-(currentHit[1]-i)])/ len(cdsSeq[currentHit[0]: i])*100
		if percentIdentity==100.0:
			hit=[currentHit[0], i, currentHit[2], currentHit[3]-(currentHit[1]-i), currentHit[4], percentIdentity]	
			return hit
		else:
			pass

	return []

def trimLeft(ALPHA_TRIM, currentHit, cdsSeq,geneSeq) :
	"""
	This function 

	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:
	ALPHA_TRIM:
	Returns
	-------
	hit :
	"""
	for i in range (currentHit[0]+ALPHA_TRIM, currentHit[1]-ALPHA_TRIM, ALPHA_TRIM ):
		
		se1=cdsSeq[i: currentHit[1]]
		se2=geneSeq[currentHit[2]+(i-currentHit[0]): currentHit[3]]
		if len(se1) >=currentHit[1]-ALPHA_TRIM:
			break
			
		percentIdentity=compute_alignment_identity(se1,se2)/ len(se1)*100
		if percentIdentity==100.0:
			hit=[i, currentHit[1], currentHit[2]+(i-currentHit[0]), currentHit[3], currentHit[4], percentIdentity]	
			return hit
		else:
			pass

	return []

def trimRightLeft(ALPHA_TRIM, currentHit, cdsSeq,geneSeq):
	"""
	This function 

	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:
	ALPHA_TRIM:
	Returns
	-------
	hit :
	"""
	k=0
	for i in range (currentHit[0]+ALPHA_TRIM, currentHit[1]-ALPHA_TRIM, ALPHA_TRIM ):
		k =+ ALPHA_TRIM 
		if i >= currentHit[1]-i :
				return []
		else:
			
			percentIdentity=compute_alignment_identity(cdsSeq[i: currentHit[1]-i],geneSeq[currentHit[2]+k: currentHit[3]-k])/ len(cdsSeq[i: currentHit[1]-i])*100
		
			if percentIdentity==100.0:
				hit=[i, currentHit[1]-i, currentHit[2]+(i-currentHit[0]), currentHit[3]-(currentHit[1]-i), currentHit[4], percentIdentity]	
				return hit
			else:
				pass
		

	return []


def extendHit(currentHit, cdsSeq,geneSeq):
	"""
	This function extend hit
	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:
	
	Returns
	-------
	extendedHit :
	"""
	hit= extendRight(currentHit, cdsSeq,geneSeq)
	
	extendedHit= extendLeft(hit, cdsSeq,geneSeq)
	
	return extendedHit
def extendRight(currentHit, cdsSeq,geneSeq):
	"""
	This function extend hit to the rigth
	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:

	Returns
	-------
	hit :
	"""
	k=0
	for i, j in zip(cdsSeq[currentHit[1]:],geneSeq[currentHit[3]:]):
	
		if i == j:
			k=k+1	
		else:
			
			break
	 
	hit= [currentHit[0], currentHit[1]+k, currentHit[2], currentHit[3]+k, currentHit[4], currentHit[5]]
	
	return hit

def extendLeft(currentHit, cdsSeq,geneSeq):
        """
	This function extend hit to the left
	Parameters
	----------
	cdsSeq:
        geneSeq:
        currentHit:

	Returns
	-------
	hit :
	"""
	k=0
	se1= reversed(cdsSeq[:currentHit[0]])
	se2= reversed(geneSeq[:currentHit[2]])
	for i, j in zip(se1,se2):
		
		if i == j:
			k=k+1	
		else:
			
			break
	hit= [currentHit[0]-k, currentHit[1], currentHit[2]-k, currentHit[3], currentHit[4], currentHit[5]]
	
	return hit
def trimAGGT(hit1,hit2, cdsSeq,geneSeq ):
	"""
	This function trim hits based on splice signals

	Parameters
	----------
	cdsSeq:
        geneSeq:
        hit1:
	hit2:
	Returns
	-------
	int
	"""
	#case one
	if geneSeq[hit1[3]:hit1[3]+2] == 'GT' and geneSeq[hit2[2]+(hit1[1]-hit2[0])-2:hit2[2]+(hit1[1]-hit2[0])]== 'AG':
		
		return 1
	#case two
	elif geneSeq[hit1[3]-(hit1[1]-hit2[0]):hit1[3]-(hit1[1]-hit2[0])+2] == 'GT' and geneSeq[hit2[2]-2:hit2[2]]== 'AG':
		
		return 2
	#other
	else:
		
		return 1

	

def extractandExtend(blocks, cdsSeq,geneSeq):
    """
    This function extract and extend remaining segment of cds
    Parameters
    ----------

    
    blocks: list
        list of blocks found by the local alignment
    cdsSeq:
    geneSeq:
    Returns
    -------
    blocks:
    remainingBlocs:
    """
    blocks=[[0,0,0,0,0,0]]+ order(   blocks) +[[len(cdsSeq), len(cdsSeq),len(geneSeq),len(geneSeq),0,0]]
    
    remainingBlocs=[]
    for i in range (0,len(blocks)-1):
	
	if blocks[i][1]==  blocks[i+1][0]:
			pass
	else:
		bloc= [blocks[i][1], blocks[i+1][0], blocks[i][3], blocks[i+1][2]]
		foundBloc= extendExon(bloc,cdsSeq,geneSeq)
		if len(foundBloc) > 0:
			blocks.append(foundBloc)
		else:
			remainingBlocs.append(bloc)
    blocks.remove([0,0,0,0,0,0])
    blocks.remove([len(cdsSeq), len(cdsSeq),len(geneSeq),len(geneSeq),0,0])
    return blocks, remainingBlocs


def extendExon(bloc,cdsSeq,geneSeq):
	"""
	This function extend segment
	Parameters
	----------


	blocks: list
	list of blocks found by the local alignment
	cdsSeq:
	geneSeq:

	Returns
	-------
	foundBloc
	"""
        foundBloc=[]
	
	sequence1= cdsSeq[bloc[0]: bloc[1]]
	sequence2= geneSeq[bloc[2]: bloc[3]]
	if len(sequence2) >= len(sequence1):
		codeSequence1= codeSegment(sequence1)
		k=0
		for j in range (bloc[2], bloc[3]-len(sequence1)+1, 1):
			if k==0:

				codeSequence2= codeSegment(geneSeq[j: j+len(sequence1)])
			else:
				codeSequence2= codeNewSegment(codeSequence2, sequence2[k-1], sequence2[k+len(sequence1)-1],  len(sequence1)-1  )
			if codeSequence1 == codeSequence2:
					foundBloc=[bloc[0], bloc[1], j,j+len(sequence1) ]
					return foundBloc
		
			k += 1
		return []
					
def codeSegment(sequence ):
	"""
	This function 

	Parameters
	----------


	sequence:
	Returns
	-------
	codeSequence:
	"""
 	dictCode= {
		"A" : 0,
		"C" : 1,
		"G" : 2,
		"T" : 3		
	}
	codeSequence=0
	for i, j in enumerate(sequence):
		
		codeSequence += dictCode[j]* (4**i)
	return codeSequence

def codeNewSegment(previousScore, previousAlphabet, currentAlphabet,  l  ):
	"""
	This function 
	Parameters
	----------

	previousScore:
        previousAlphabet:
        currentAlphabet:
        l :
	Returns
	-------
	codeSequence:
	"""
 	dictCode= {
		"A" : 0,
		"C" : 1,
		"G" : 2,
		"T" : 3		
	}
	codeSequence=0
	
	codeSequence = ((previousScore - dictCode[previousAlphabet]) /4 )+(dictCode[currentAlphabet]*(4**l))
	return codeSequence

def launch_tblastx(cdsfile, genefile,evalue, block1,block2,cdsid, geneid):
	"""
	This function 
	Parameters
	----------

	cdsfile:
	genefile:
	evalue:
	block1:
	block2:
	cdsid:
	geneid:

	Returns
	-------
	blastoutput:
	"""
	block1_qs,block1_qe,block1_ss,block1_se=block1
	block2_qs,block2_qe,block2_ss,block2_se=block2

	blastoutput = str(os.getcwd()+'/src/results/blast_results/CDS_'+ cdsid+'_'+str(block1_qe)+'-'+str(block2_qs)+'_VS_Gene_'+ geneid+'_'+str(block1_se)+'-'+str(block2_ss)+'.tblastx')
    
	command = "tblastx -query " + cdsfile + " -subject " + genefile +  " -evalue " + str(evalue) + " -outfmt "+ str(7) +  " -strand plus " +  " -out " + str(blastoutput)
#              " qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq  -strand plus "

	'''
	tblastx[-h][-help][-import_search_strategy
	filename]
	[-export_search_strategy filename][-db
	database_name]
	[-dbsize num_letters][-gilist
	filename] [-seqidlist filename]
	[-negative_gilist filename][-entrez_query
	entrez_query]
	[-db_soft_mask filtering_algorithm][-db_hard_mask
	filtering_algorithm]
	[-subject subject_input_file][-subject_loc
	range] [-query input_file]
	[-out output_file][-evalue
	evalue] [-word_size int_value]
	[-qcov_hsp_perc float_value][-max_hsps
	int_value]
	[-xdrop_ungap float_value][-searchsp
	int_value] [-sum_stats bool_value]
	[-max_intron_length length][-seg
	SEG_options]
	[-soft_masking soft_masking][-matrix
	matrix_name]
	[-threshold float_value][-culling_limit
	int_value]
	[-best_hit_overhang float_value][-best_hit_score_edge
	float_value]
	[-window_size int_value][-lcase_masking][-query_loc
	range]
	[-strand strand][-parse_deflines][-query_gencode
	int_value]
	[-db_gencode int_value][-outfmt
	format] [-show_gis]
	[-num_descriptions int_value][-num_alignments
	int_value]
	[-line_length line_length][-html][-max_target_seqs
	num_sequences]
	[-num_threads int_value][-remote][-version]
	'''
	os.system(command)

	return blastoutput

def concatenate_exons(blocklist,cdsseq,geneseq):
	"""
	This function concatenates all pairs of successive blocks that 
	whose cds_separation is less than MAX_CDS_SEPARATION_FOR_CONCAT, and 
	gene_separation doesnot exeed MAX_CDS_SEPARATION_FOR_CONCAT +/-MAX_SEPARATION_DIFF_FOR_CONCAT.
	Parameters
	----------

	blocklist:
	cdsseq:
	geneseq:

	Returns
	-------
	blocklist:
	"""


	block_to_delete=[]

	blocklist = order(blocklist)
    
	for i in range(len(blocklist)-1):

		block1_qs,block1_qe,block1_ss,block1_se=blocklist[i]
		block2_qs,block2_qe,block2_ss,block2_se=blocklist[i+1]

		cdsseparation = block2_qs-block1_qe
		geneseparation = block2_ss-block1_se

		if (cdsseparation <= MAX_CDS_SEPARATION_FOR_CONCAT):
		    
		    diff = abs(cdsseparation - geneseparation)
		    if (diff <= MAX_SEPARATION_DIFF_FOR_CONCAT and diff%3 == 0):
		        
		        alignment = pairwise2.align.globalms(geneseq[block1_se-5:block2_ss+5], cdsseq[block1_qe-5:block2_qs+5],2.0, 0.0, -5.0, -1.0)
		        alignementscore = -1 
		        for aln in alignment:
		            alignementscore = aln[4] 
		            
		        if (alignementscore > 0):
		            blocklist[i+1]=[block1_qs,block2_qe,block1_ss,block2_se]
		            block_to_delete.append(blocklist[i])

	for block in block_to_delete:
		blocklist.remove(block)

    	return blocklist

