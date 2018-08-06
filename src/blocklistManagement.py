#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``blocklistmanagement.py`` **module description**:

This modules contains all the function used to manage blocklists.

.. moduleauthor:: Safa Jammali and Jean-David Aguilar

2017-2018

"""

from ete3 import *
from Bio import pairwise2

CANONICAL_DONOR = "GT"
CANONICAL_ACCEPTOR = "AG"

NC_DONOR1="GT"
NC_ACCEPTOR1="AG"
NC_DONOR2="GT"
NC_ACCEPTOR2="AG"


def cover_wholeCDS(blocklist,cdslength):
    """
    This function returns True if the blocklist cover the entire CDS,
    and False otherwise.
    
    Parameters
    ----------

    blocklist:
    cdslength:

    Returns
    -------
    boolean:
    """
    percentcover = cover_percentage(blocklist,cdslength)
    if (percentcover < 100.0):
        return False
    else:
        return True

def is_continuous_onCDS(blocklist):
	"""
	This function 

	Parameters
	----------

	blocklist:


	Returns
	-------
	boolean:
	"""
	for i in range(1,len(blocklist)):
		if(blocklist[i-1][1] < blocklist[i][0]):
		    return False
	return True

def cover_percentage(blocklist,cdslength):
    """
    This function returns the percentage of coverage of the cds 
    by the blocklist.

    Parameters
    ----------

    blocklist:
    cdslength:

    Returns
    -------
    percentcover:
    """
    percentcover = 0.0
    
    coveredlength = 0

    for i in range(len(blocklist)):
        block = blocklist[i] 
        coveredlength += block[1]-block[0]
        if(i > 0):
            prev_block =  blocklist[i-1]
            if(block[2]-prev_block[3] == 0):
                coveredlength += block[0]-prev_block[1]
        
    percentcover= 100.0 * coveredlength / cdslength
    return percentcover

def order(blocklist):
    """
    This function orders blocks in the blocklist by increasing
    order of query start location.
    Parameters
    ----------

    blocklist:
    

    Returns
    -------
    blocklist:
   
    """

    for i in range(len(blocklist)):
        min = i
        for j in range(i+1,len(blocklist)):

            if(blocklist[j][0] < blocklist[min][0]):
                min = j
        if(min != i):
            tmp = blocklist[min]
            blocklist[min] = blocklist[i]
            blocklist[i] = tmp
    return blocklist


def launch_nwpairwise(sequence1, sequence2):
    """
    This function launches global pairwise alignment
    Parameters
    ----------

    sequence1:
    sequence2:
    

    Returns
    -------
    alignment:
   
    """
    alignment = pairwise2.align.globalms(sequence1, sequence2,2,0,-5,-1)
    return alignment

def compute_block_identity(cds, gene,block):
    """
    This function is use to launch a global alignment between the 2 current sequence 
    and to supply the compute_alignment_identity function with the 2 aligned sequence 
    Parameters
    ----------

    cds:
    gene:
    block:
    

    Returns
    -------
    block_identity:
    """
    block_identity = 0.0
    block_qs = block[0] #query start
    block_qe = block[1] #query end
    block_ss = block[2] #subject start
    block_se = block[3] #subject end

    cds_= cds[block_qs:block_qe]
    gene_= gene[block_ss:block_se]
    sequence1 = ""
    sequence2 = ""

    if(len(cds_)!=0 and len(gene_)!=0):
        
        if(len(cds_)==len(gene_)):
            sequence1 = gene_
            sequence2 = cds_
        else:
            # maximise matches and minimize gaps
            alignment = launch_nwpairwise(gene_, cds_)
            sequence1, sequence2 = alignment[0][0],alignment[0][1]
        
        block_identity = 1.0 * compute_alignment_identity(sequence1, sequence2) /len(sequence2)

    return block_identity

def compute_alignment_identity(sequence1, sequence2):
    """
    This function aim to return the percentage of identity of 2 sequences
     Parameters
    ----------

    sequence1:
    sequence2:
    

    Returns
    -------
    aln_identity:
    
    """

    aln_identity = 0.0
    match = 0
    length = len(sequence1)
   
    for i in range(length):
        if(sequence1[i] == sequence2[i]):
            match += 1
    aln_identity = 1.0 * match
    return aln_identity


def create_target_cds(blocklist,geneseq):
    """
    This function 
    Parameters
    ----------

    blocklist:
    geneseq:
    

    Returns
    -------
    targetcds:
    
    """
    targetcds = ""
    
    for i in range(len(blocklist)):
        block = blocklist[i]
        targetcds +=  geneseq[block[2]:block[3]]

    return targetcds


def create_target_exon(blocklist):
    """
    This function 
    Parameters
    ----------

    blocklist:
        

    Returns
    -------
    targetexon:
    
    """
    targetexon = []

    for i in range(len(blocklist)):
        block = blocklist[i]
        targetexon.append(block[2:])

    return targetexon

def is_protein(targetcds):
    """
    This function return True if the predicted transcript have these 3 protein feature : length multiple of 3, start and stop codon
    Parameters
    ----------

    targetcds:
        

    Returns
    -------
    boolean:
    
    """

    length = len(targetcds)
        
    if(length%3 == 0 and targetcds[:3] == 'ATG' and (targetcds[length-3:] == 'TAA' or targetcds[length-3:] == 'TAG' or targetcds[length-3:] == 'TGA')):
        return True
    else:
        return False
    
def compute_splicing_sites(blocklist,geneseq):
	"""
	This function 
	Parameters
	----------

	blocklist:
	geneseq:


	Returns
	-------
	splicing_sites:

	"""
	splicing_sites = []
	for i in range(len(blocklist)-1):
		block1 = blocklist[i]
		block2 = blocklist[i+1]
		splicing_sites.append([geneseq[block1[3]:block1[3]+2],
			      geneseq[block2[2]-2:block2[2]]])
	return splicing_sites

def is_valid_splicing(splice):
    """
    This function 
    Parameters
    ----------

    
    splice:
        

    Returns
    -------
    
    
    """
    donor, acceptor = splice
    return ((donor == CANONICAL_DONOR and acceptor == CANONICAL_ACCEPTOR) or
            (donor == NC_DONOR1 and acceptor == NC_ACCEPTOR1) or
            (donor == NC_DONOR2 and acceptor == NC_ACCEPTOR2))
    
def all_valid_splicing_sites(splicing_sites):
    """
    This function 
    Parameters
    ----------

    
    splicing_sites:
        

    Returns
    -------
    all-valid
    
    """
    all_valid = True
    for splice in splicing_sites:
        all_valid = all_valid and is_valid_splicing(splice)
    return all_valid

def print_blocklist(cdsid, geneid, cds, gene, blocklist,outfile, outputformat):
    """
    This function is the main function of the graphic 
    representation of the prédiction.
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    blocklist:
    outfile: 
    outputformat:
        

    Returns
    -------
   
    
    """
    cds_len = len(cds)

    if(len(blocklist) > 0):
        first_block = blocklist[0]
        first_qs = first_block[0]
        if(0 < first_qs):
            string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [0,first_qs], outputformat)
            outfile.write(string_to_print)

        
        for i in range(len(blocklist)-1):
            block = blocklist[i]
            string_to_print = compute_aln_string(cdsid, geneid, cds, gene, block, outputformat)
            outfile.write(string_to_print)

            next_block = blocklist[i+1]
            block_qe = block[1]
            next_block_qs = next_block[0]
            if(block_qe < next_block_qs):
                string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [block_qe ,next_block_qs], outputformat)
                outfile.write(string_to_print)

        last_block = blocklist[-1]
        string_to_print = compute_aln_string(cdsid, geneid, cds, gene, last_block, outputformat)
        outfile.write(string_to_print)

        last_block_qe = last_block[1]
        if(last_block_qe < cds_len):
            string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [last_block_qe,cds_len], outputformat)
            outfile.write(string_to_print)
        
    else:
        string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [0 ,cds_len], outputformat)
        outfile.write(string_to_print)
        
    outfile.write("\n")
        

def print_wholealignment(cdsid, geneid, cds, gene, blocklist,outfile, outputformat):
    """
    This function
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    blocklist:
    outfile: 
    outputformat:
        

    Returns
    -------
   
    
    """
    cds_len = len(cds)
    gene_len = len(cds)
    cds_=""
    gene_ = ""

    sequence1 = ""
    sequence2 = ""
    if(len(cds_)==len(gene_)):
        sequence1 = gene_
        sequence2 = cds_
    elif(len(cds_)== 0):
        sequence1 = gene_
        sequence2 = '-' * len(sequence1)
    elif(len(gene_)== 0):
        sequence2 = cds_
        sequence1 = '-' * len(sequence2)
    else:
        alignment = pairwise2.align.globalms(gene_, cds_,2,0,-10,-1)
        sequence1, sequence2 = alignment[0][0],alignment[0][1]

        
    for i in range(len(blocklist)):
        block = blocklist[i]

        if (i==0):# first block
            cds_ += '-' * block[2]
            gene_ += gene[:block[2]]
            cds_ += cds[:block[0]]            
            gene_ += '-' * block[0]

        else:
            cds_ += '-' * (block[2] - blocklist[i-1][3])
            gene_ += gene[blocklist[i-1][3]:block[2]]
            cds_ += cds[blocklist[i-1][1]:block[0]]            
            gene_ += '-' * (block[0] - blocklist[i-1][1])
            
        sequence1 = ""
        sequence2 = ""
        cdsblock = cds[block[0]:block[1]]
        geneblock = gene[block[2]:block[3]]
        if(len(cdsblock)==len(geneblock)):
            sequence1 = geneblock
            sequence2 = cdsblock
        elif(len(cdsblock)== 0):
            sequence1 = geneblock
            sequence2 = '-' * len(sequence1)
        elif (len(geneblock)==0) :
            sequence2 = cdsblock
            sequence1 = '-' * len(sequence2)
        else:
            
            alignment = pairwise2.align.globalms(geneblock, cdsblock,2,0,-10,-1)
            sequence1, sequence2 = alignment[0][0],alignment[0][1]
        gene_ += sequence1
        cds_ += sequence2

        cds_ += '-' * (gene_len - blocklist[-1][3])
        gene_ += gene[blocklist[-1][3]:gene_len]
        cds_ += cds[blocklist[-1][1]:cds_len]            
        gene_ += '-' * (cds_len - blocklist[-1][1])

    outfile.write(">"+geneid+"\n")
    outfile.write(gene_+"\n")
    outfile.write(">"+cdsid+"\n")
    outfile.write(cds_+"\n")

          
def compute_aln_string(cdsid, geneid, cds, gene,block, outputformat):
    """
    This function produce the visual representation for each aligned block using a global alignment
     
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    block:
     
    outputformat:
        

    Returns
    -------
    string_to_print
    
    """
    string_to_print = ""
    
    block_qs = block[0] #query start
    block_qe = block[1] #query start
    block_ss = block[2] #subject start
    block_se = block[3] #subject end
    block_identity = "%.2f" % (compute_block_identity(cds, gene,block))
    gene_= gene[block_ss:block_se]
    cds_= cds[block_qs:block_qe]

    sequence1 = ""
    sequence2 = ""
    if(len(cds_)==len(gene_)):
        sequence1 = gene_
        sequence2 = cds_
    elif(len(cds_)== 0):
        sequence1 = gene_
        sequence2 = '-' * len(sequence1)
    elif(len(gene_)== 0):
        sequence2 = cds_
        sequence1 = '-' * len(sequence2)
    else:
        alignment = pairwise2.align.globalms(gene_, cds_,2,0,-10,-1)
        sequence1, sequence2 = alignment[0][0],alignment[0][1]

    aln_length = len(sequence1)

    string_to_print = cdsid + "\t" + geneid + "\t" + str(aln_length) + "\t" + str(block_qs) + "\t" + str(block_qe) + "\t" + str(block_ss) +  "\t" + str(block_se) + "\t" + str(block_identity) +  "\t" + gene[block_ss-2:block_ss] + "<Exon>" + gene[block_se:block_se+2] + "\n"
    
    if(outputformat == "aln"):
        sequence1 = gene[block_ss-BORDER_LENGTH:block_ss] + sequence1 + gene[block_se:block_se+BORDER_LENGTH]
        sequence2 = BORDER_LENGTH*" " + sequence2 + BORDER_LENGTH*" "

        aln_srspair = format_alignment(sequence1,sequence2)

        string_to_print +=  aln_srspair
        
    return string_to_print
  
def compute_unaln_string(cdsid, geneid, cds, gene,interval, outputformat):
    """
    This function produce the visual representation for each unaligned cds block
       
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    interval:
     
    outputformat:
        

    Returns
    -------
    string_to_print
    
   
    """
    string_to_print = ""
    
    block_qs = interval[0] #query start
    block_qe = interval[1] #query end
    cds_= cds[block_qs:block_qe]

    length = block_qe - block_qs
    
    string_to_print = cdsid + "\t" + geneid + "\t" + str(length) + "\t" + str(block_qs) + "\t" + str(block_qe) + "\t" + "-" +  "\t" + "-" + "\n"
    
    if(outputformat == "aln"):

        sequence = format_sequence(cds_)

        string_to_print +=  sequence
        
    return string_to_print

def format_alignment(sequence1,sequence2):
    """
    This function is use to produce the visualization format with terminal adapt line length
       
    Parameters
    ----------
    sequence1:
    sequence2:
     
           

    Returns
    -------
    aln_srspair:
    
    
    """

    aln_srspair = ""

    markup_line = compute_markup_line(sequence1,sequence2)
    
    start = 0
    while start < len(sequence1):
        line_length = LENGTH_ALN_LINE
        if(start+line_length > len(sequence1)):
            line_length = len(sequence1)-start
        subsequence1 = sequence1[start:start+line_length]
        subsequence2 = sequence2[start:start+line_length]
        submarkup_line = markup_line[start:start+line_length]
        
        aln_srspair +=  subsequence1 + "\n"
        aln_srspair +=  submarkup_line + "\n"
        aln_srspair += subsequence2  + "\n\n"
        start += line_length

    return aln_srspair

def format_sequence(seq):
    """
    This function 

    Parameters
    ----------
    seq:
     
           

    Returns
    -------
    sequence:
    
    
    """
    sequence = ""
    
    start = 0
    while start < len(seq):
        line_length = LENGTH_ALN_LINE
        if(start+line_length > len(seq)):
            line_length = len(seq)-start
        subseq = seq[start:start+line_length]
        
        sequence +=  subseq + "\n"
        start += line_length

    sequence += "\n"
    return sequence


def compute_markup_line(sequence1,sequence2):
    """
    This function is use to return the central line of the visual alignment with the match and mismatch
    Parameters
    ----------
    sequence1:
    sequence2:
         

    Returns
    -------
        
    markup_line:
    
    """

    markup_line = ""
    for i in range(len(sequence1)):
        if(sequence1[i]==sequence2[i]):
            markup_line += "|"
        else:
            markup_line += " "
    return markup_line

def print_exonextremitylist(cdsid, geneid, blocklist, cdsexonendlistfile,geneexonstartlistfile,geneexonendlistfile):
    """
    This function 

    Parameters
    ----------
    cdsid:
    geneid:
    blocklist:
    cdsexonendlistfile:
    geneexonstartlistfile:
    geneexonendlistfile:
         

    Returns
    -------
        
    
    """
    cdsexonendlistfile.write(">"+cdsid+"\n")
    geneexonstartlistfile.write(">"+geneid+"\n")
    geneexonendlistfile.write(">"+geneid+"\n")
    for block in  blocklist:
        cdsexonendlistfile.write(str(block[1])+" ")
        geneexonstartlistfile.write(str(block[2])+" ")
        geneexonendlistfile.write(str(block[3])+" ")
    cdsexonendlistfile.write("\n")
    geneexonstartlistfile.write("\n")
    geneexonendlistfile.write("\n")
