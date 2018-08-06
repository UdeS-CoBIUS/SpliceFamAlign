#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``Exact.py`` **module description**:

This module is the module that launches the comparison for all pairs of source CDS and target gene using the SFA_G method.

.. moduleauthor:: Safa Jammali



2017-2018

"""
import blocklistManagement
from blocklistManagement import *
import alignment
from alignment import *
import numpy
from numpy import *


LIMIT_ADD_NUCLEOTIDES = 3
STATUS_EXISTING_PROTEIN = 1
STATUS_PREDICTED_PROTEIN = 2
STATUS_COMPLETE_CDS = 3
STATUS_PARTIAL_CDS = 4

MIN_IDENTITY_FINAL = 0.5
ThresholdIdentityMatch = 30
GLOBALLENGTH= 100


GAP= -700
MATCH=  1000
MISMATCH= -1044

RealIntron = 1

minimumIntron =50
maximumIntron=5000
CORRECT_EXON =20


BLAST_ERROR_ALIGN=5

INTERVAL_BLOC_COMPATIBLE = 5
INTERVAL_BLOC_COMPATIBLE_GENE_CDS = 40
REAL_SPLICE_SITES =6000
ONE_REAL_SPLICE_SITES = 3000

REAL_EXON_JUNCTION = 6000

SPLICE_SIGNAL_GTAG = -4270
SPLICE_SIGNAL_GCAG = -5314
SPLICE_SIGNAL_ATAC = -6358
SPLICE_SIGNAL_OTHER = -7395


MIN_EVALUE = 1.0/10000000 # 10-7
#MIN_EVALUE = 1.0/100 # 10-2
###############################
### LOCAL COMPARISON ##########
###############################

def localAlignment(cdsId,cdsLen, cdsExon,geneId, geneLen, geneSeq):
	"""
	This function calls launchAndTrimTblastx function to do local comparison between CDS and gene and returns blocks, that are
	  significant detected hits between CDS and gene

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


    Returns
      -------
      blockList: list
          bloks and information about the alignment of each pairs of CDS against gene
      """





	blockList = []

	# give e value
	evalue = MIN_EVALUE

	# fins cds file's name and gene file's name
	cdsFile = os.getcwd() + '/src/sequences/cds/CDS_' + cdsId + '.fasta'
	geneFile = os.getcwd() + '/src/sequences/genes/Gene_' + geneId + '.fasta'

	# calls tblastx to do alignment
	blockList = launchAndTrimTblastx(cdsFile, geneFile, evalue, cdsId, geneId,cdsLen, geneLen, cdsExon, geneSeq)
	if (len(blockList) == 0):
		print "No alignment was found with tblastx"
		return []

	# returns list of hits if found, and returns [] if not found
	return blockList

def launchAndTrimTblastx( cdsFile, geneFile, evalue, cdsId, geneId,cdsLen, geneLen, cdsExon, geneSeq):
	"""
	This function  extracIntervalBlock block to align, call localAlignmentTblastx to do the local alignment and returns blocks, that are
	  significant detected hits between CDS and gene

	Parameters
	----------

	geneId: string
		gene id
	cdsId: string
		cds id
	cdsExon: list
		list of exons of cds
	cdsFile: file name
		cds file name that contains the cds sequence
	geneFile: file name
		gene file name that contains the gene sequence
	evalue: int
		cutt off value  for blast significant result
	geneLen:int
		length of gene
	cdsLen:int
		length of cds

    Returns
	-------
	blockList: list
		bloks and information about the alignment of each pairs of CDS against gene
	"""


	blastOutput = launch_tblastx(cdsFile, geneFile, evalue, [0, 0, 0, 0],[cdsLen, cdsLen, geneLen, geneLen], cdsId, geneId)

	blockList = parseTblastxTrimBlock(cdsId, geneId, cdsExon, blastOutput, 0, 0, geneSeq)
	os.remove(blastOutput)
	blockList = order(blockList)

	# returns list of hits
	return blockList

def parseTblastxTrimBlock(cdsId, geneId, cdsExon, blastOutput, blockQueryEnd, blockSubjectEnd, geneSeq):

	"""
	This function takes tblastx results hits, treats compatibility and returns blocks that are compatible

	Parameters
	----------

	blockQueryEnd : int
			value of the end residu to adjust blast result on the querry
	blockSubjectEnd: int
			value of the end residu to adjust blast result on the subject
	geneId: string
		gene id
	cdsId: string
		cds id
	blastOutput: filename
		filename that contains tblastx results (hits informations)
	cdsExon: list
		list of exons of cds

	Returns
	-------
	blocks: list
		bloks and information about the alignment of each pairs of CDS against gene
	"""

	
	dictExonBlocks = {}
	dictExonBlocksToaffect = {}
	incompatibledictExonBlocksToaffect = {}
	compatibledictExonBlocksToaffect = {}
	blocks = []
	valueExonBlocksNew = []

	
	firstExon= cdsExon[0]
	lastExon = cdsExon[len(cdsExon)-1]

	listExonOfCDS =[]
	#create an empty dictionary for each exon of cds
	for bloc in cdsExon:
		listExonOfCDS.append([bloc[0],bloc[1]])
		dictExonBlocks[(bloc[0],bloc[1])] =  []
		dictExonBlocksToaffect[(bloc[0],bloc[1])] =  []
		incompatibledictExonBlocksToaffect[(bloc[0],bloc[1])] =  []
		compatibledictExonBlocksToaffect[(bloc[0], bloc[1])] = []
	listHits=[]
	'''
	with open(blastOutput, 'rb') as csvfile:
		spamreader = csv.reader(csvfile,delimiter = ' ', quotechar = ' ')

		for row in spamreader:


			if row[0][0] is not '#':
				# construct a hit


				beginNewHitCDS= int(row[0].split('\t')[6])+int(blockQueryEnd) - 1
				endNewHitCDS = int(row[0].split('\t')[7]) + int(blockQueryEnd)
				beginNewHitgene = int(row[0].split('\t')[8])+ int(blockSubjectEnd) - 1
				endNewHitgene = int(row[0].split('\t')[9]) + int(blockSubjectEnd)
				evalueNewHit = float(row[0].split('\t')[10])
				if (beginNewHitCDS < endNewHitCDS and beginNewHitgene < endNewHitgene and evalueNewHit <= MIN_EVALUE):
					listHits.append([beginNewHitCDS,endNewHitCDS, beginNewHitgene, endNewHitgene, evalueNewHit])
	'''
        listHits= readBlastOut(blastOutput, blockQueryEnd, blockSubjectEnd)

	listHits = order(listHits)
	
	for hit in listHits:
				# construct a hit

				beginNewHitCDS= hit[0]
				endNewHitCDS = hit[1]
				beginNewHitgene = hit[2]
				endNewHitgene = hit[3]
				evalueOfNewHit = hit[4]
	
				#affect hit in dictExonBlocks

				for indiceExonOfCDS in range (0, len(listExonOfCDS)):
					res1 = beginNewHitCDS >= listExonOfCDS[indiceExonOfCDS][0]
					res2 = beginNewHitCDS < listExonOfCDS[indiceExonOfCDS][1]
					res= res1 and res2

					if res:
						listExonOfCDS = listExonOfCDS[indiceExonOfCDS:len(listExonOfCDS)]


						break
					else:
						pass

	
				#  determine exon that contains the new hit
				maxLengthIntersection = -1
				segmentInteraction= None
				beginExonChoice = -1
				endExonChoice = -1
				i = 0
				continuer = True
				
				while (continuer == True and i <= len(listExonOfCDS)-1):
					e = listExonOfCDS[i]

					[beginExon, endExon] = e
					if beginExon >= endNewHitCDS:
						continuer == False
						break



					segmentInteractionLast, lengthIntersection = determineIntersection(beginExon, endExon,\
						int(beginNewHitCDS), int(endNewHitCDS), int(beginNewHitgene), int(endNewHitgene),cdsId, geneId,\
																					   geneSeq, firstExon, lastExon)


					if lengthIntersection > maxLengthIntersection:
							maxLengthIntersection = lengthIntersection
							segmentInteraction = segmentInteractionLast

							[beginExonChoice , endExonChoice]= e
					i = i+1
				
				dictExonBlocks[(beginExonChoice, endExonChoice)].append((segmentInteraction, evalueOfNewHit))
	
	for exontoTreat in cdsExon:
		
		[beginExon, endExon] = exontoTreat
		listHitsExon = dictExonBlocks[(beginExon, endExon)]
		
		if len(listHitsExon) == 0 :
			pass
		elif len(listHitsExon) == 1 :
			
			[(NewHit, evalueOfNewHit)] = listHitsExon

			dictExonBlocksToaffect[(beginExon, endExon)].append((NewHit, evalueOfNewHit))
		else:

			listHitsExonInOrder =sorted(listHitsExon, key=lambda colonnes: colonnes[1])

			
			for hitExon in listHitsExonInOrder:
				

				(NewHit, evalueOfNewHit) = hitExon
				if len(dictExonBlocksToaffect[(beginExon, endExon)]) == 0:
					dictExonBlocksToaffect[(beginExon, endExon)].append((NewHit, evalueOfNewHit))

				else:
					
					dictExonBlocksToaffect[(beginExon, endExon)] = affectHit(dictExonBlocksToaffect[(beginExon, endExon)], NewHit, cdsId, geneId, evalueOfNewHit)

	
	for i in range(0, len(dictExonBlocksToaffect.keys())):
		
		for j in range(i + 1, len(dictExonBlocksToaffect.keys())):
			
			for hitexon in dictExonBlocksToaffect.values()[i]:
				exon1 = dictExonBlocksToaffect.keys()[i]
				for hitotherexon in dictExonBlocksToaffect.values()[j]:
					exon2 = dictExonBlocksToaffect.keys()[j]
					
					(hit1, eval1) = hitexon
					(hit2, eval2) = hitotherexon
					if exon1[1] <= exon2[0]:
						

						if int(hit1[3]) > int(hit2[2]):
							
							if (min(float(eval1), float(eval2)) == eval1):
								
								incompatibledictExonBlocksToaffect[exon2].append(hitotherexon)  #
							else:
								
								incompatibledictExonBlocksToaffect[exon1].append(hitexon)
					if exon2[1] <= exon1[0]:
						

						if int(hit2[3]) > int(hit1[2]):
							
							if (min(float(eval1), float(eval2)) == eval1):
								incompatibledictExonBlocksToaffect[exon2].append(
									hitotherexon)  
							else:
								incompatibledictExonBlocksToaffect[exon1].append(
									hitexon) 
	for cle, val in incompatibledictExonBlocksToaffect.items():
		
		if len(val) == 0:
			
			if len(dictExonBlocksToaffect[cle]) != 0:
				for k in dictExonBlocksToaffect[cle]:
					compatibledictExonBlocksToaffect[cle].append(k)
		
		else:
			newVal = []
			for j in dictExonBlocksToaffect[cle]:
				newVal = dictExonBlocksToaffect[cle]

			for k in val:
				if k in newVal:
					newVal.remove(k)
			compatibledictExonBlocksToaffect[cle] = newVal
		

	for blockItem in compatibledictExonBlocksToaffect.values():
			
			if (len(blockItem) == 0):
				pass

			if(len (blockItem) == 1):


				[(hitToAdd, eval) ]= blockItem
				blocks.append(hitToAdd)
			else:

				for element in blockItem:

					(hitToAdd, eval) = element
					blocks.append(hitToAdd)
	
	return blocks



def determineIntersection(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, firstExon, lastExon):
	"""
	This function determine intersection between cds exon and hit, it returns the intersection with its length

	Parameters
	----------

	beginExon: int
		value of begin exon
	endExon: int
		value of end exon
	beginNewHit: int
		value of begin of hit in cds
	endNewHit: int
		value of end of hit in cds
	beginNewHitgene: int
		value of begin of hit in gene
	endNewHitgene: int
		value of end of hit in gene

	geneId: string
		gene id
	cdsId: string
		cds id


	Returns
	-------

	newhit: list
		list that contains the new information of intersection between hit and exon
	length: int
		value of the length of the intersection between hit and exon

	"""

	length = 0
	newhit = None

	if beginNewHit > endExon :
		pass

	elif beginExon > endNewHit:
		pass

	#if newHit in exon
	elif (beginExon <= beginNewHit and beginNewHit <= endExon) and (beginExon <= endNewHit and endNewHit <= endExon):


		verifDonor = verifyLengthAndSitesDonor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, lastExon)
		if ((endExon - endNewHit ) <= 	LIMIT_ADD_NUCLEOTIDES) or verifDonor:
			endInCDS = endExon
			endInGene = endNewHitgene + (endExon - endNewHit)

		else:
			endInCDS = endNewHit
			endInGene = endNewHitgene

		verifAcceptor = verifyLengthAndSitesAcceptor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, firstExon)
		if (abs(beginNewHit - beginExon) <= LIMIT_ADD_NUCLEOTIDES) or verifAcceptor:
			beginInCDS = beginExon
			beginInGene =  beginNewHitgene - (beginNewHit - beginExon)

		else:
			beginInCDS = beginNewHit
			beginInGene = beginNewHitgene

		newhit = [beginInCDS, endInCDS , beginInGene, endInGene]
		length = endInCDS - beginInCDS

		
	# if exon in newHit
	elif beginNewHit <= beginExon and beginExon <= endNewHit and beginNewHit <= endExon and endExon <= endNewHit:

		beginNewHitGeneLevel = beginNewHitgene + abs(beginExon - beginNewHit)
		endNewHitGeneLevel = endNewHitgene - abs(endNewHit - endExon)
		newhit = [beginExon, endExon, beginNewHitGeneLevel,endNewHitGeneLevel]
		length = endExon - beginExon


		
	# if beginExon <= beginNewHit <= endExon  < endNewHit
	elif ((beginExon <= beginNewHit and beginNewHit <= endExon) and endNewHit > endExon):

		endNewHitGeneLevel = endNewHitgene - abs(endNewHit - endExon)
		verifAcceptor = verifyLengthAndSitesAcceptor(beginExon, endExon, beginNewHit, endExon, beginNewHitgene,
													 endNewHitGeneLevel, cdsId, geneId, geneSeq, firstExon)
		if (abs(beginNewHit - beginExon) <= LIMIT_ADD_NUCLEOTIDES) or verifAcceptor :
			beginInCDS = beginExon
			beginInGene = beginNewHitgene - abs(beginNewHit - beginExon)

		else:
			beginInCDS = beginNewHit
			beginInGene = beginNewHitgene

		newhit = [beginInCDS, endExon, beginInGene, endNewHitGeneLevel]
		length = endExon - beginInCDS
		
	# if beginNewHit < beginExon  <= endNewHit  <= endExon
	elif (beginNewHit < beginExon and (beginExon <= endNewHit and endNewHit <= endExon)):
		beginNewHitGeneLevel = beginNewHitgene + abs(beginExon - beginNewHit)

		verifDonor = verifyLengthAndSitesDonor(beginExon, endExon, beginExon, endNewHit, beginNewHitGeneLevel,
											   endNewHitgene, cdsId, geneId, geneSeq, lastExon)


		if ((endExon - endNewHit) <= LIMIT_ADD_NUCLEOTIDES) or verifDonor:

			endInCDS= endExon
			endNewHitGeneLevel = endNewHitgene + abs(endExon - endNewHit)

		else:
			endInCDS = endNewHit
			endNewHitGeneLevel = endNewHitgene

		newhit = [beginExon, endInCDS, beginNewHitGeneLevel, endNewHitGeneLevel]
		length = endInCDS - beginExon
		

	else:
		
		print 'Verif case didn take account in  this function determineIntersection'\
			, beginExon,endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene
		exit(-1)
	return newhit, length


def verifyLengthAndSitesAcceptor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, firstExon):
	"""
    This function verify the length of hit , if it's > than exon length and the acceptor is AG , it returns True if not it returns false

    Parameters
    ----------

    beginExon: int
        value of begin exon
    endExon: int
        value of end exon
    beginNewHit: int
        value of begin of hit in cds
    endNewHit: int
        value of end of hit in cds
    beginNewHitgene: int
        value of begin of hit in gene
    endNewHitgene: int
        value of end of hit in gene

    geneId: string
        gene id
    cdsId: string
        cds id

	geneSeq: string
    	gene sequence
    Returns
    -------

    boolean

    """

	lenExon =  endExon - beginExon
	lenHit = endNewHit - beginNewHit
	hitCoverExon = False
	if lenHit >= (lenExon/2):
		hitCoverExon = True



	if beginExon > beginNewHit:

		endAcceptor = beginNewHitgene + abs(beginNewHit - beginExon)
	elif beginNewHit > beginExon:
		endAcceptor = beginNewHitgene - abs(beginNewHit - beginExon)

	else:
		endAcceptor = beginNewHitgene
	beginAcceptor = endAcceptor -2

	acceptor = geneSeq[beginAcceptor: endAcceptor] #acceptor== 'AG' end intron
	if  (((hitCoverExon and acceptor == 'AG') == True) or ((hitCoverExon and ([beginExon, endExon] == firstExon)) == True)):
		return True
	else:
		return False



def verifyLengthAndSitesDonor(beginExon, endExon, beginNewHit, endNewHit, beginNewHitgene, endNewHitgene, cdsId, geneId, geneSeq, lastExon):
	"""
    This function verify th length of hit if it s > than exon length and the donor is GT, it returns True if not it returns false

    Parameters
    ----------

    beginExon: int
        value of begin exon
    endExon: int
        value of end exon
    beginNewHit: int
        value of begin of hit in cds
    endNewHit: int
        value of end of hit in cds
    beginNewHitgene: int
        value of begin of hit in gene
    endNewHitgene: int
        value of end of hit in gene

    geneId: string
        gene id
    cdsId: string
        cds id

	geneSeq: string
    	gene sequence
    Returns
    -------

    boolean

    """

	lenExon =  endExon - beginExon
	lenHit = endNewHit - beginNewHit

	hitCoverExon = False
	if lenHit >= (lenExon/2):
		hitCoverExon = True



	if endExon > endNewHit:
		beginDonor = abs(endExon - endNewHit) +endNewHitgene

	elif endNewHit > endExon:
		beginDonor = endNewHitgene - abs(endExon - endNewHit)
	else:
		beginDonor = endNewHitgene

	endDonor = beginDonor + 2


	donor = geneSeq[beginDonor: endDonor]  #donor== 'GT' start intron

	if ( ((hitCoverExon and donor == 'GT') == True )or ((hitCoverExon and ([beginExon, endExon] == lastExon)) == True)):
		return True
	else:
		return False


def affectHit(dictExonBlocksNew, segmentInteraction, cdsId, geneId, evalueOfNewHit):
	"""
	This function  compares and treats compatibility and returns  the compatible hits

	Parameters
	----------

	dictExonBlocksNew: list
		list of compatibl hits
	segmentInteraction:list
		new hit that will be treat

	geneId: string
		gene id
	cdsId: string
		cds id
	evalueOfNewHit: int
		value of e value given by tblastx for the hit

	Returns
	-------

	dictExonBlocksNew: list
		 contains the update list of hits


	"""

	dictExonBlocksNew_ = []
	countHit = 0
	HitAndEvalue = 0
	
	while HitAndEvalue < len(dictExonBlocksNew):
		
		(currentHit, evalueHit) = dictExonBlocksNew[HitAndEvalue]

		
		minBeginCds = min(segmentInteraction[0], currentHit[0])
		minBeginGene = min(segmentInteraction[2], currentHit[2])

		maxBeginCds = max(segmentInteraction[0], currentHit[0])
		maxBeginGene = max(segmentInteraction[2], currentHit[2])

		maxEndCds = max(segmentInteraction[1], currentHit[1])
		maxEndGene = max(segmentInteraction[3], currentHit[3])

		minDistanceCDS = abs(minBeginCds - maxBeginCds)
		minDistanceGene = abs(minBeginGene - maxBeginGene)
		compatible = 2
		case = ''
		# currentHit is include in the newHit, verify compatibility and add newHit or not
		if segmentInteraction[0] <= currentHit[0] and currentHit[0] <= segmentInteraction[1] \
				and segmentInteraction[0] <= currentHit[1] and currentHit[1] <= segmentInteraction[1] \
				and segmentInteraction[2] <= currentHit[2] and currentHit[2] <= segmentInteraction[3] \
				and segmentInteraction[2] <= currentHit[3] and currentHit[3] <= segmentInteraction[3]:
			
			countHit = countHit + 1
			# verify compatibility and add newHit or not
			if abs(minDistanceCDS - minDistanceGene) <= INTERVAL_BLOC_COMPATIBLE_GENE_CDS:

				# add the newHit: segmentInteraction
				if countHit == len(dictExonBlocksNew):
					
					dictExonBlocksNew_.append((segmentInteraction, evalueOfNewHit))
					
					break
				else:
				
					HitAndEvalue = HitAndEvalue + 1
				
			else:
				
				if (min(evalueHit, evalueOfNewHit) == evalueHit):
					dictExonBlocksNew_.append((currentHit, evalueHit))
					break

				else:
					
					if countHit == len(dictExonBlocksNew):

						
						dictExonBlocksNew_.append((segmentInteraction, evalueOfNewHit))
						break
					# remove currentHit
					else:

						HitAndEvalue = HitAndEvalue + 1
					
		# currentHit is different to the newHit, add newHit
		elif (segmentInteraction[0] > currentHit[1] and segmentInteraction[2] > currentHit[3]) \
				or (currentHit[0] > segmentInteraction[1] and currentHit[2] > segmentInteraction[3]):
			

			countHit = countHit + 1
			if countHit == len(dictExonBlocksNew):
				
				dictExonBlocksNew_.append((segmentInteraction, evalueOfNewHit))
				dictExonBlocksNew_.append((currentHit, evalueHit))
				break
			else:
				dictExonBlocksNew_.append((currentHit, evalueHit))
				HitAndEvalue = HitAndEvalue + 1
			

		# currentHit and newHit are overlaping, verify compatibility and add newHit or not
		elif (minBeginCds == segmentInteraction[0] and minBeginGene == segmentInteraction[2]) \
				or (minBeginCds == currentHit[0] and minBeginGene == currentHit[2]):
			
			countHit = countHit + 1
			compatible = 1
		else:
			countHit = countHit + 1
			
			compatible = 0
			case = 'non compatibe'

		if compatible == 1:

			if abs(minDistanceCDS - minDistanceGene) <= INTERVAL_BLOC_COMPATIBLE_GENE_CDS:

				
				if countHit == len(dictExonBlocksNew):
					
					dictExonBlocksNew_.append(
						([minBeginCds, maxEndCds, minBeginGene, maxEndGene], min(evalueHit, evalueOfNewHit)))
					break
				else:
					
					segmentInteraction = [minBeginCds, maxEndCds, minBeginGene, maxEndGene]
					evalueOfNewHit = min(evalueHit, evalueOfNewHit)
					HitAndEvalue = HitAndEvalue + 1

				
			else:
				
				if (min(evalueHit, evalueOfNewHit) == evalueHit):
					dictExonBlocksNew_.append((currentHit, evalueHit))
					break
				else:
					

					if countHit == len(dictExonBlocksNew):
						
						dictExonBlocksNew_.append((segmentInteraction, evalueOfNewHit))
						break
					else:

						HitAndEvalue = HitAndEvalue + 1

		elif (compatible == 0 and case == 'non compatibe'):

			if (min(evalueHit, evalueOfNewHit) == evalueHit):
				dictExonBlocksNew_.append((currentHit, evalueHit))
				break
			
			else:
				
				if countHit == len(dictExonBlocksNew):
					
					dictExonBlocksNew_.append((segmentInteraction, evalueOfNewHit))
					break
				else:

					HitAndEvalue = HitAndEvalue + 1

				
	return dictExonBlocksNew_

#########################
### alignement PD ##########
#########################


def InterBlocksPD(resultsBloc,  blockListDP, geneExon,cdsExon, intronList,  cdsId, geneId, lenGene,lenCDS, cdsSeq, geneSeq ):
	"""
	This function extracts the inter blocks and launchs exact comparison for each interblock using dynamic prog.
	At the end, it recovers blocks and interblocks in the same list

	Parameters
	----------

	resultsBloc: list
			list contains blocks ound by local alignment eof gene and CDS sequences
	geneExon: list
        list of exon of gene:
	cdsExon: list
		list of exons of cds
	intronList:dictionary
		dictionary of gene intron, the keys are gene Id and the values are the blocks interval of intron gene
	cdsId: string
		CDS Id
	geneId: string
		gene Id
	lenGene: int
		length of gene sequence
	lenCDS: int
		length of CDS sequence
	cdsSeq: string
		   CDS sequence
	lenGene:  string
			gene sequence


	Returns
	-------
	allBlocks: list
		list of bloks found for each pairs of CDS against gene
	"""

	
	newBlocs=[]

	exonInGene ={}
	intronInGene ={}
	exonBlocksWithinCds ={}
	posExonStartGene=[]
	posExonEndGene=[]

	posExonStartCds=[]
	posExonEndCds=[]
	posStartIntron=[]
	posEndIntron =[]

	
	extendedBlocks= order(blockListDP)
	# for each blocks
	
	for bloci in extendedBlocks:


		#identify extremity of two successifs blokcs
		[QueryStart, QueryEnd, SubjectStart, SubjectEnd]= bloci
		
		
		if (QueryEnd -QueryStart)* (SubjectEnd - SubjectStart)> GLOBALLENGTH:
			
			pass
		else:

			
			if  ( QueryEnd -QueryStart ) <= 0 or(SubjectEnd - SubjectStart ) <=0 :
				print   cdsId, ' ', geneId, ' len(segmentCDS) or len(segmentgene) == 0 \n '

			elif (QueryEnd -QueryStart ) > 0 and (SubjectEnd - SubjectStart ) >0:

				

				exonInGene, intronInGene, exonBlocksWithinCds,posExonStartGene, posExonEndGene, posExonStartCds, posExonEndCds, posStartIntron, posEndIntron \
					= cdsGenerestrictionExonIntron(geneId, cdsId, geneExon, cdsExon, intronList, QueryStart, QueryEnd, \
												   SubjectStart, SubjectEnd)

				
				newBloc = DynamicProgrammationAlignment(geneId, cdsId, geneSeq, cdsSeq,
														QueryStart, QueryEnd, \
														SubjectStart, SubjectEnd,\
														geneExon, cdsExon, posExonStartGene,
														posExonEndGene, posExonStartCds,
														posExonEndCds, posStartIntron, posEndIntron, exonInGene, intronInGene, exonBlocksWithinCds)

				
				#add the newu block found
				for  bloc in newBloc:

					newBlocs += [bloc]





	# take all blocks
	allBlocks = resultsBloc + newBlocs

	#order blocks
	allBlocks= order(allBlocks)
	

	allBlocks = concatenate_exons(allBlocks, cdsSeq, geneSeq)
	
	return allBlocks


def cdsGenerestrictionExonIntron(geneId, cdsId,  geneExon, cdsExon, intronList, cdsBegin, cdsEnd, geneBegin, geneEnd):
	"""
	This function extracts the blocks of exons and introns of gene and exons of CDS belonging to the inter block

	Parameters
	----------

	geneId: string
		gene Id
	cdsId: string
		cds Id

	geneExon: list
		list of exon of gene
	cdsExon: list
		list of exons of cds

	intronList:dictionary
			   dictionary of gene intron, the keys are gene Id and the values are the blocks interval of intron gene

	cdsBegin: int
		value of the begining of the interblock at cds level
	cdsEnd: int
		end value  of the interblock at cds level
	geneBegin: int
		value of the beginig of the interblock at gene level
	geneEnd:int
		end value  of the interblock at gene level
	Returns
	-------

	posNucleotideInExonStart: dictionary
		dictionary of start positions of nucleotide in exon of gene belonging the interblock
	posNucleotideInExonEnd: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock
	posNucleotideInExonStartCds: dictionary
		dictionary of start positions of nucleotide in exon of cds  belonging the interblock
	posNucleotideInExonEndCds: dictionary
		dictionary of end positions of nucleotide in exon of cds  belonging the interblock
	posStartIntron: dictionary
		dictionary of start positions of nucleotide in intron  belonging the interblock
	posEndIntron:dictionary
		dictionary of end positions of nucleotide in intron  belonging the interblock

	exonInGene: dictionary
		 dictionary of list of exons in each gene, taking account the limit of interblock
	intronInGene: dictionary
		dictionary of list of introns in each gene, taking account the limit of interblock
	exonBlocksWithinCds: dictionary
		dictionary of list of exons in each cds, taking account the limit of interblock



	"""

	exonInGene={}
	posNucleotideInExonStart={}
	posNucleotideInExonEnd={}

	intronInGene = {}
	posStartIntron = {}
	posEndIntron = {}

	exonBlocksWithinCds = {}
	posNucleotideInExonStartCds = {}
	posNucleotideInExonEndCds = {}

	# extract interval block of exons  of gene belonging to the inter block
	exonInGene, posNucleotideInExonStart, posNucleotideInExonEnd =extracIntervalBlock(geneExon, geneId, exonInGene, \
												posNucleotideInExonStart, posNucleotideInExonEnd, geneBegin, geneEnd)

	# extract interval block of introns  of gene  belonging to the inter block
	intronInGene, posStartIntron , posEndIntron =extracIntervalBlock(intronList, geneId, intronInGene, posStartIntron,\
																	 posEndIntron, geneBegin, geneEnd)

	# extract interval block of exons  of CDS belonging to the inter block
	exonBlocksWithinCds, posNucleotideInExonStartCds, posNucleotideInExonEndCds, =extracIntervalBlock(cdsExon, cdsId, \
					exonBlocksWithinCds, posNucleotideInExonStartCds, posNucleotideInExonEndCds, cdsBegin, cdsEnd)


	#return
	return  exonInGene, intronInGene, exonBlocksWithinCds, posNucleotideInExonStart, posNucleotideInExonEnd, \
			posNucleotideInExonStartCds, posNucleotideInExonEndCds, posStartIntron, posEndIntron




def extracIntervalBlock(sourceDictionary, sourceId, blockInInterval, posNucleotideStartInInterval,\
						posNucleotideEndInInterval, beginInterval, endInterval):
	"""
	This function computes the blocks of exons and introns of gene and exons of CDS belonging to the inter block

	Parameters
	----------

	sourceDictionary: dictionary
		dictionary that contains list of cds exons or gene exons or gene introns
		 the keys are source Id and the values are the blocks interval of exon/intron source
	sourceId: string
		it can be cds id or gene Id

	blockInInterval: dictionary
		dictionary that contains list of cds exons or gene exons or gene introns beloning to the inter block
	posNucleotideStartInInterval: dictionary
		dictionary that cotains list of begin postion of cds exons or gene exons or gene introns beloning to the inter block
	posNucleotideEndInInterval: dictionary
		dictionary that cotains list of end postion of cds exons or gene exons or gene introns beloning to the inter block
	beginInterval: int
		begin of the inter block
	endInterval: int
		end of the inter block

	Returns
	-------
	blockInInterval:  dictionary
		 dictionary of list of exons/introns
	posNucleotideStartInInterval: dictionary
		dictionary of start positions of nucleotide
	posNucleotideEndInInterval: dictionary
		dictionary of end positions of nucleotide


	"""
	
	for key, values in sourceDictionary.items():

		if sourceId == key:

			blockInInterval[sourceId] = []
			posNucleotideStartInInterval[sourceId] = []
			posNucleotideEndInInterval[sourceId] = []
			for i in range(0, len(values)):
				minInter = values[i][0]
				maxInter = values[i][1]
				posNucleotideStartInInterval[sourceId] += [minInter]
				posNucleotideEndInInterval[sourceId] += [maxInter]

				if beginInterval > maxInter:
					pass
				if minInter > endInterval:
					pass
				if beginInterval <= minInter and minInter <= endInterval and beginInterval <= maxInter and maxInter <= endInterval:
					blockInInterval[sourceId].append([minInter, maxInter])
				if beginInterval <= minInter and minInter < endInterval and maxInter > endInterval:
					blockInInterval[sourceId].append([minInter, endInterval])
				if beginInterval > minInter and beginInterval < maxInter and maxInter <= endInterval:
					blockInInterval[sourceId].append([beginInterval, maxInter])

				if beginInterval > minInter and maxInter > endInterval:
					blockInInterval[sourceId].append([beginInterval, endInterval])





	return blockInInterval, posNucleotideStartInInterval, posNucleotideEndInInterval

def DynamicProgrammationAlignment(geneId, cdsId, geneSequence, cdsSequence, \
								  cdsBeginSegment, cdsEndSegment, geneBeginSegment, geneEndSegment, \
								  geneExon, cdsExon, posExonStartGene, posExonEndGene, \
								  posExonStartCds, posExonEndCds, posStartIntron, posEndIntron, exonInGene, intronInGene,\
								  exonBlocksWithinCds):





	"""
	This function do a dynamic programmation alignment. It has  cases to compute and to choose the maximum one
	to fill it in the DP table

	Parameters
	----------

	cdsId: string
		CDS Id
	geneId: string
		gene Id
	cdsSequence: string
		   CDS sequence
	geneSequence:  string
		   gene sequence
	geneExon: dictionary
		dictionay with gene id as key and exons interval list as value
	cdsExon: dictionary
		dictionay with CDS id as key and exons interval list as value
	cdsBeginSegment: int
		value of begin segment in cds
	cdsEndSegment: int
		value of end segment in cds
	geneBeginSegment: int
		value of begin segment in gene
	geneEndSegment: int
		value of end segment in gene

	posExonStartGene: dictionary
		dictionary of start positions of nucleotide in exon of gene belonging the interblock
	posExonEndGene: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock
	posExonStartCds: dictionary
		dictionary of start positions of nucleotide in exon of cds  belonging the interblock
	posExonEndCds: dictionary
		dictionary of end positions of nucleotide in exon of cds  belonging the interblock
	posStartIntron: dictionary
		dictionary of start positions of nucleotide in intron  belonging the interblock
	posEndIntron:dictionary
		dictionary of end positions of nucleotide in intron  belonging the interblock

	exonInGene: dictionary
		 dictionary of list of exons in each gene, taking account the limit of interblock
	intronInGene: dictionary
		dictionary of list of introns in each gene, taking account the limit of interblock
	exonBlocksWithinCds: dictionary
		dictionary of list of exons in each cds, taking account the limit of interblock

	Returns
	-------
	newBlocs: list
		list of blocks exact aligned
	"""



	segmentCDS = '-' + cdsSequence[cdsBeginSegment: cdsEndSegment]
	segmentgene = '-' + geneSequence[geneBeginSegment: geneEndSegment]

	# create DP table
	M = create_matrix(len(segmentCDS), len(segmentgene))
	# create trace table
	Mtrace = create_matrix(len(segmentCDS), len(segmentgene))


	# initilization
	for i in range(1, len(segmentgene)):
		M[i][0] = 0
		Mtrace[i][0] = [0, 0, 'casei']

	for j in range(1, len(segmentCDS)):
		M[0][j] = j * GAP
		Mtrace[0][j] = [0, 0, 'casej']

	# filling the DP table
	for j in range(1, len(segmentCDS)):
		for i in range(1, len(segmentgene)):

			i1 = 0
			j1 = 0
			iPrime = 0
			iMinusL = 0

			# calcul six cases in DP
			case1, case2, case3,  case5, iPrime, case6, iMinusL = case_alignment(geneSequence, \
																										  cdsSequence,
																										  M, i, j,
																										  geneExon,
																										  cdsExon,
																										  geneId, cdsId,
																										  posExonStartGene,
																										  posExonEndGene,
																										  posExonStartCds, \
																										  posExonEndCds,
																										  cdsBeginSegment,
																										  geneBeginSegment,
																										  segmentCDS,
																										  segmentgene, \
																										  posStartIntron,
																										  posEndIntron,
																										   exonInGene, \
																				intronInGene, exonBlocksWithinCds)




			# choose maximum case of the six
			M[i][j] = max(case1, case2, case3, case5, case6)
			# fill in the box according to the chosen case and store its trace
			if M[i][j] == case1:
				Mtrace[i][j] = [i - 1, j - 1, 'case1']
			elif M[i][j] == case2:
				Mtrace[i][j] = [i , j-1, 'case2']
			elif M[i][j] == case3:
				Mtrace[i][j] = [i-1, j , 'case3']
			elif M[i][j] == case5:
				Mtrace[i][j] = [iPrime, j, 'case5']
			elif M[i][j] == case6:
				Mtrace[i][j] = [iMinusL, j, 'case6']

	posmaxligne = foundPosition(M,segmentCDS, segmentgene )

	# display, identification and knowledge of the blocks alignment constructed by PD
	newBloc = traceAlignmentDP(posmaxligne, len(segmentCDS) - 1, Mtrace, geneId, cdsId, cdsBeginSegment,
							   geneBeginSegment, geneSequence, cdsSequence)

	# returns aligned blocks
	return newBloc


def create_matrix(m, n):
	"""


	Parameters
	----------

	m: int
		lengtn of cds
	n: int
		length of gene
	Returns
	-------
	M: list
		matrix n*m

	"""
	M = [[0 for i in range(m)] for j in range(n)]
	return M


def case_alignment(geneSequence, cdsSequence, M, i, j, geneExon, cdsExon, geneId, cdsId, posExonStartGene, posExonEndGene, \
				   posStartExonOnCds, posEndExonOnCds, cdsBeginSegment, geneBeginSegment, sequenceCdsInterblocks,\
				   sequenceGeneInterblocks, posStartIntron, posEndIntron, exonInGene, intronInGene, exonBlocksWithinCds):

	"""
	This function computes cases given in DP algorithms.

	Parameters
	----------

	cdsId: string
		CDS Id
	geneId: string
		gene Id
	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence

	geneExon: dictionary
		dictionay with gene id as key and exons interval list as value
	cdsExon: dictionary
		dictionay with CDS id as key and exons interval list as value

	cdsBeginSegment: int
		value of begin segment in cds
	geneBeginSegment: int
		value of begin segment in gene
	sequenceCdsInterblocks: string
		    CDS sequence in the interblock
	sequenceGeneInterblocks: string
		   origin CDS sequence in the interblock

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	posExonStartGene: dictionary
		dictionary of start positions of nucleotide in exon of gene belonging the interblock
	posExonEndGene: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock
	posExonStartCds: dictionary
		dictionary of start positions of nucleotide in exon of cds  belonging the interblock
	posExonEndCds: dictionary
		dictionary of end positions of nucleotide in exon of cds  belonging the interblock
	posStartIntron: dictionary
		dictionary of start positions of nucleotide in intron  belonging the interblock
	posEndIntron:dictionary
		dictionary of end positions of nucleotide in intron  belonging the interblock

	exonInGene: dictionary
		 dictionary of list of exons in each gene, taking account the limit of interblock
	intronInGene: dictionary
		dictionary of list of introns in each gene, taking account the limit of interblock
	exonBlocksWithinCds: dictionary
		dictionary of list of exons in each cds, taking account the limit of interblock

	Returns
	-------
	case1: int
		score of case1
	case2:int
		score of case2
	case3:int
		score of case3
	maxcas5: int
		score of case5
	iPrime: int
		position in gene, beging of real intron
	maxcas6:int
		score of case6
	iMinusL: int
		position in gene, beging of detected intron
	"""

	
	verifRealExonJunction = verifRealExonJunctionOnCDS(j+cdsBeginSegment-1, posStartExonOnCds[cdsId])


	
	# Compute case 5 score
	if  len(intronInGene[geneId]) == 0 :

		maxcas5 = -Infinity
		iPrime = 0

	else:
		maxcas5, iPrime = maxIntronscore(M, i, j, verifRealExonJunction, intronInGene[geneId], geneSequence, cdsSequence,\
						sequenceGeneInterblocks, sequenceCdsInterblocks, cdsBeginSegment, geneBeginSegment,\
						posExonStartGene[geneId], posExonEndGene[geneId], posStartIntron[geneId], posEndIntron[geneId])


	
	# Compute case 6 score
	maxcas6, iMinusL = maxIntronscore_case6(M, i, j, verifRealExonJunction, sequenceGeneInterblocks, sequenceCdsInterblocks,\
											geneBeginSegment, posStartIntron[geneId], posEndIntron[geneId])


	
	# Compute case 1 score
	case1 = case_1(M, i, j, MATCH, MISMATCH, sequenceGeneInterblocks, sequenceCdsInterblocks)
	# Compute case 2 score
	case2 = case_2(M, i, j, GAP)
	# Compute case 3 score
	case3 = case_3(M, i, j, GAP)
	
	#return scores
	return case1, case2, case3,  maxcas5, iPrime, maxcas6, iMinusL


def verifRealExonJunctionOnCDS(j,posStartExonOnCds ):

	"""
	This boolean function verif if actual postion in cds is a real junction or not.

	Parameters
	----------

	j: int
		actual position in CDS

	posExonStartCds: list
		list of start positions of nucleotide in exon of cds  belonging the interblock

	Returns
	-------
	boolean : True if yes ano False if not

	"""
	if j in posStartExonOnCds:
		 return True
	return False

def maxIntronscore(M, i, j, verifRealExonJunction, listIntronGene, geneSequence, cdsSequence, sequenceGeneInterblocks, \
				   sequenceCdsInterblocks, cdsBeginSegment, geneBeginSegment, posExonStartGene, posExonEndGene,\
				   posStartIntron, posEndIntron):


	"""
	This function computes score of case 5

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	verifRealExonJunction: boolean
		True if actual postion in cds is a real junction

	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence

	cdsBeginSegment: int
		value of begin segment in cds
	geneBeginSegment: int
		value of begin segment in gene
	sequenceCdsInterblocks: string
		    CDS sequence in the interblock
	sequenceGeneInterblocks: string
		   origin CDS sequence in the interblock

	posExonStartGene: dictionary
		dictionary of start positions of nucleotide in exon of gene  belonging the interblock

	posExonEndGene: dictionary
		dictionary of end positions of nucleotide in exon of gene  belonging the interblock

	posStartIntron: list
		list of start positions of nucleotide in exon of gene  belonging the interblock
	posEndIntron: list
		list of end positions of nucleotide in exon of gene  belonging the interblock


	listIntronGene: list
		 list of introns in each gene, taking account the limit of interblock

	Returns
	-------
	maxCsoreCase5: int
		score of case5
	iPrime: int
		position in gene, beging of real intron

	"""

	maxScoreCase5= -Infinity
	iPrime = 0

	scoreCase5 = -Infinity

	# compute score of real exon junction
	if verifRealExonJunction:
		scoreRealExonJunction = REAL_EXON_JUNCTION

	else:
		scoreRealExonJunction = 0

	# compute score of case 5 and its gene intron begining
	scoreSpliceSignals = SPLICE_SIGNAL_OTHER
	for intron in range(0, len(listIntronGene)):
		startIntron = listIntronGene[intron][0]
		endIntron = listIntronGene[intron][1]
		if i + geneBeginSegment -1 == endIntron  :

			scoreSpliceSites = verifRealSpliceSites(startIntron, endIntron, posStartIntron, posEndIntron)

			donor = sequenceGeneInterblocks[startIntron - geneBeginSegment : startIntron - geneBeginSegment + 2]
			acceptor = sequenceGeneInterblocks[endIntron - geneBeginSegment - 2 : endIntron - geneBeginSegment]

			if donor== 'GT' and  acceptor == 'AG':
				scoreSpliceSignals = SPLICE_SIGNAL_GTAG
			elif donor == 'GC' and acceptor == 'AG':
				scoreSpliceSignals = SPLICE_SIGNAL_GCAG
			elif donor == 'AT' and acceptor == 'AC':
				scoreSpliceSignals = SPLICE_SIGNAL_ATAC

			scoreCase5 = M[startIntron - geneBeginSegment ][j] + scoreRealExonJunction + scoreSpliceSignals + scoreSpliceSites

			if maxScoreCase5 < scoreCase5:
				maxScoreCase5 = scoreCase5
				iPrime = startIntron - geneBeginSegment +1


		else:
			pass
		
	#return
	return maxScoreCase5, iPrime



def verifRealSpliceSites(i, l, posStartIntron, posEndIntron):

	"""
	This function verify if positions in gene are real start and real end intron in the gene
	and based on thes , it returns a score

	Parameters
	----------


	i: int
		it can be a start intron
	l: int
		it can be end of intron

	posStartIntron: list
		list of start positions of nucleotide in exon of gene  belonging the interblock
	posEndIntron: list
		list of end positions of nucleotide in exon of gene  belonging the interblock


	Returns
	-------
	scoreRealSpliceSites : int
		score of splice sites, it can be 0, 1000 or 2000

	"""
	scoreRealSpliceSites = 0
	if i in posStartIntron:
		scoreRealSpliceSites += ONE_REAL_SPLICE_SITES

	if l in  posEndIntron:
		scoreRealSpliceSites += ONE_REAL_SPLICE_SITES

	#return scores
	return scoreRealSpliceSites


def maxIntronscore_case6(M, i, j, verifRealExonJunction, genesequence, CDSsequence, geneBeginSegment, posStartIntron, posEndIntron):

	"""
	This function computes score of case 6

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	verifRealExonJunction: boolean
		True if actual postion in cds is a real junction

	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence
	geneBeginSegment: int
		value of begin segment in gene



	posStartIntron: list
		list of start positions of nucleotide in exon of gene  belonging the interblock
	posEndIntron: list
		list of end positions of nucleotide in exon of gene  belonging the interblock


	Returns
	-------
	maxcas6: int
		score of case 6
	iMinusL: int
		position in gene, begining of intron


	"""

	maxcas6 = -Infinity
	iMinusL = 0
	scoreCase6 = -Infinity

	# compute score of real exon junction
	if verifRealExonJunction:
		scoreRealExonJunction = REAL_EXON_JUNCTION

	else:
		scoreRealExonJunction = 0

	# compute score of case 6 and its gene intron begining
	scoreSpliceSignals = SPLICE_SIGNAL_OTHER

	scoreSpliceSites = -Infinity
	donor = ''
	acceptor = ''
	if i + geneBeginSegment > minimumIntron and i + geneBeginSegment < maximumIntron:

		if i > minimumIntron:
			for l in range(i - minimumIntron, 1, -1):

				
				donor = genesequence[l + geneBeginSegment - 1: l + geneBeginSegment + 1]
				acceptor = genesequence[i + geneBeginSegment - 1: i + geneBeginSegment + 1]

				if donor == 'GT' and acceptor == 'AG':
					scoreSpliceSignals = SPLICE_SIGNAL_GTAG
				elif donor == 'GC' and acceptor == 'AG':
					scoreSpliceSignals = SPLICE_SIGNAL_GCAG
				elif donor == 'AT' and acceptor == 'AC':
					scoreSpliceSignals = SPLICE_SIGNAL_ATAC

				scoreSpliceSites = verifRealSpliceSites(i + geneBeginSegment - 1, l + geneBeginSegment - 1,
														posStartIntron, posEndIntron)

				scoreCase6 = M[l][j] + scoreRealExonJunction + scoreSpliceSignals + scoreSpliceSites

				if maxcas6 < scoreCase6:
					maxcas6 = scoreCase6
					iMinusL = l

			

	elif i + geneBeginSegment > maximumIntron:
		
		if i > minimumIntron and i < maximumIntron:
			for l in range(i - minimumIntron, 1, -1):

				
				donor = genesequence[l + geneBeginSegment - 1: l + geneBeginSegment + 1]
				acceptor = genesequence[i + geneBeginSegment - 1: i + geneBeginSegment + 1]

				if donor == 'GT' and acceptor == 'AG':
					scoreSpliceSignals = SPLICE_SIGNAL_GTAG
				elif donor == 'GC' and acceptor == 'AG':
					scoreSpliceSignals = SPLICE_SIGNAL_GCAG
				elif donor == 'AT' and acceptor == 'AC':
					scoreSpliceSignals = SPLICE_SIGNAL_ATAC

				scoreSpliceSites = verifRealSpliceSites(i + geneBeginSegment - 1, l + geneBeginSegment - 1,
														posStartIntron, posEndIntron)

				scoreCase6 = M[l][j] + scoreRealExonJunction + scoreSpliceSignals + scoreSpliceSites

				if maxcas6 < scoreCase6:
					maxcas6 = scoreCase6
					iMinusL = l
				
		elif i > maximumIntron:
			for l in range(i - minimumIntron, i - maximumIntron, -1):

				
				donor = genesequence[l + geneBeginSegment - 1: l + geneBeginSegment + 1]
				acceptor = genesequence[i + geneBeginSegment - 1: i + geneBeginSegment + 1]

				if donor == 'GT' and acceptor == 'AG':
					scoreSpliceSignals = SPLICE_SIGNAL_GTAG
				elif donor == 'GC' and acceptor == 'AG':
					scoreSpliceSignals = SPLICE_SIGNAL_GCAG
				elif donor == 'AT' and acceptor == 'AC':
					scoreSpliceSignals = SPLICE_SIGNAL_ATAC

				scoreSpliceSites = verifRealSpliceSites(i + geneBeginSegment - 1, l + geneBeginSegment - 1,
														posStartIntron, posEndIntron)

				scoreCase6 = M[l][j] + scoreRealExonJunction + scoreSpliceSignals + scoreSpliceSites

				if maxcas6 < scoreCase6:
					maxcas6 = scoreCase6
					iMinusL = l
			


	else:
		maxcas6 = -Infinity
		iMinusL = 0
		scoreCase6 = -Infinity

	#return
	return maxcas6, iMinusL



def case_1(M, i, j,  MATCH, MISMATCH, genesequence, CDSsequence):
	"""
	This function computes score of case 1

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	MATCH: int
		score of matching residu
	MISMATCH: int
		score of matching residu

	cdsSequence: string
		   origin CDS sequence
	geneSequence:  string
		   origin gene sequence

	Returns
	-------
	res: int
		score of case 2
	"""


	if genesequence[i] == CDSsequence[j]:
		score = MATCH
	else:
		score = MISMATCH
	res = M[i - 1][j - 1] + score

	return res


def case_2(M, i, j, GAP):
	"""
	This function computes score of case 2

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	GAP: int
		score of gap


	Returns
	-------
	res: int
		score of case 2



	"""
	res = M[i ][j- 1] + GAP

	return res



def case_3(M, i, j, GAP):
	"""
	This function computes score of case 2

	Parameters
	----------

	M: list
		DP matrix
	i: int
		position in gene
	j: int
		position in cds

	GAP: int
		score of gap


	Returns
	-------
	res: int
		score of case 3



	"""
	res = M[i- 1][j ] + GAP

	return res



def foundPosition(M,segmentCDS, segmentgene):
	"""
	This function computes the position that contains the max in the DP matrix M

	Parameters
	----------


	M: list
		DP matrix

	segmentCDS: string
		CDS sequence
	segmentgene: string
		gene sequence



	Returns
	-------
	posmaxligne: int
		position that contains the max
	"""

	# found max coord
	fincds = []
	



	for i in range(0, len(segmentgene)):

			fincds.append(M[i][len(segmentCDS) - 1])


	'''
	for j in range(0, len(segmentCDS)):

		for i in range(0, len(segmentgene)):

			if j == len(segmentCDS) - 1:
				fincds.append(M[i][j])
			#if i == len(segmentgene) - 1:
				#fingene.append(M[i][j])
	'''

	posmaxligne = 0

	

	maxfinc = -Infinity


	for j in range(0, len(fincds)):
		if fincds[j] > maxfinc:
			maxfinc = fincds[j]
	
	return posmaxligne

def traceAlignmentDP(positionLine, positionColumn, M, geneId, cdsId, cdsBeginSegment, geneBeginSegment,geneSequence, cdsSequence):
	"""
	This function computes the list of blocks by doing the trace back on the DP matrix

	Parameters
	----------

	cdsId: string
		CDS Id
	geneId: string
		gene Id
	M: list
		DP matrix

	positionLine:

	positionColumn:


	cdsBeginSegment: int
		value of begin segment in cds

	geneBeginSegment: int
		value of begin segment in gene


	Returns
	-------
	newBlocs: list
		list of blocks exact aligned
	"""

	newBloc = []
	casenewBloc = []



	ancien =[]
	
	while (positionLine > 0 or positionColumn > 0):

		next =  M[positionLine][positionColumn]

		if ancien == next:
			print ' infinity loop with geneId  ', geneId, '   and cdsId ', cdsId, '   ', ancien , '\n'
			break
		else:
			ancien = M[positionLine][positionColumn]

		[w, q, case] = M[positionLine][positionColumn]

		if w == 0 and q == 0 and case == 'casei':

			casenewBloc += [[positionLine, positionColumn, w, q, case]]
			break
		if w == 0 and q == 0 and case == 'casej':


			casenewBloc += [[positionLine, positionColumn, w, q, case]]
			break
		if w == positionLine - 1 and q == positionColumn - 1 and case == 'case1':


			casenewBloc += [[positionLine, positionColumn, w, q, case]]
			positionLine = positionLine - 1
			positionColumn = positionColumn - 1

		elif w == positionLine and q == positionColumn - 1  and case == 'case2':


			casenewBloc += [[positionLine, positionColumn, w, q, case]]


			positionColumn = positionColumn - 1
		elif w == positionLine- 1  and q == positionColumn  and case == 'case3':


			casenewBloc += [[positionLine, positionColumn, w, q, case]]

			positionLine = positionLine - 1

		

		elif case == 'case5':

			casenewBloc += [[positionLine, positionColumn, w, q, case]]

			positionLine = w - 1

		elif case == 'case6':



			casenewBloc += [[positionLine, positionColumn, w, q, case]]

			positionLine = w - 1

		else:
			print "-------------echec", geneId, cdsId, positionLine, positionColumn, M[positionLine][positionColumn]
			exit(-1)

		if positionLine <= 0 or positionColumn <= 0:
			casenewBloc += [[positionLine, positionColumn, w, q, case]]

			break
	
	precedent = 'stop'

	for blo in range(0, len(casenewBloc)):

		if casenewBloc[blo][4] != 'case5' and casenewBloc[blo][4] != 'case6' and casenewBloc[blo][1] != 1 and \
						casenewBloc[blo][1] != 0:
			if precedent == 'stop':
				finGene = casenewBloc[blo][0] +geneBeginSegment
				finCDS = casenewBloc[blo][1]+ cdsBeginSegment
				precedent = 'continue'

			else:
				precedent = 'continue'

		else:
			if precedent == 'stop':

				pass
			else:
				debutGene = casenewBloc[blo][0]+geneBeginSegment
				debutCDS = casenewBloc[blo][1]+ cdsBeginSegment
				newBloc+= [[debutCDS, finCDS, debutGene, finGene]]
				precedent = 'stop'

	newBlocAccepted = []

	for blocToTest in newBloc:
		[beginInCDS, endInCDS, beginInGene, endInGene] = blocToTest
		if beginInCDS != 0:
			beginInCDS = beginInCDS - 1
		if beginInGene != 0:
			beginInGene = beginInGene - 1
		if endInCDS == len(cdsSequence):
			endInCDS = endInCDS - 1
		if endInGene == len(geneSequence):
			endInGene = endInGene - 1

		sequence1 = geneSequence[beginInGene: endInGene]
		sequence2 = cdsSequence[beginInCDS: endInCDS]
		
		if len(sequence2)<=3:
			pass
		else:

			alignment = launch_nwpairwise(sequence1, sequence2)
			
			identityMatch = computeAlignmentPercentIdentity(alignment[0][0], alignment[0][1])
			
			if identityMatch >= ThresholdIdentityMatch:
				newBlocAccepted += [blocToTest]
	# returns list of blocks

	return newBlocAccepted



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
