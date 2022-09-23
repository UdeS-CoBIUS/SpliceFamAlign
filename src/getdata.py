#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``getdata.py`` **module description**:

This module is module that parses the input arguments.

.. moduleauthor:: Aïda Ouangraoua, Safa Jammali and Jean-David Aguilar

2017-2018

"""

from Bio import SeqIO
from cogent.db.ensembl import HostAccount
from cogent.db.ensembl import Compara
from cogent.db.ensembl import Species
from cogent.db.ensembl import Genome
import cogent.db.ensembl as ensembl
from blocklistManagement import *
ENSEMBL_VERSION = '85'
EXTENSION = 6000


#########################
## DATA ACQUISITION #####
#########################

def get_cds_data(gene):
    cds_data = []

    geneloc = gene.Location

    genestart = geneloc.Start - EXTENSION  
    geneend = geneloc.End + EXTENSION
    strand=  geneloc.Strand
    geneLen= geneend-genestart
    
    for transcript in gene.Transcripts:
        seq = str(transcript.Cds)
        if (transcript.BioType == 'protein_coding' and is_protein(seq)):
            new = True
            for cds in cds_data:
                if(cds[1] == seq):
                    new = False
            if(new):
                exonlist = []
                start = 0
                for i in range(len(transcript.TranslatedExons)):
                    exon = transcript.TranslatedExons[i]
                    location = str(exon.Location).split(":")[3].split("-")

                    end = start + int(location[1]) - int(location[0])
                    gstart = int(location[0]) - int(genestart)
                    gend = int(location[1]) - int(genestart)
                    
                    if strand == -1:
                        ngstart =   geneLen- gend
                        ngend =  geneLen- gstart
                       
                        exonlist.append([int(start), int(end),  int(ngstart), int(ngend)])

                        start = end
                    else:
                        ngstart = gstart
                        ngend = gend
                       
                        exonlist.append([int(start), int(end),  int(ngstart), int(ngend)])
                        start = end
               
                cds_data.append([transcript.StableId, seq, gene.StableId, exonlist])

    
    return cds_data

def get_gene_data(gene):
    gene_data = []
    gene_data.append([gene.StableId, str(getSequence(gene))])
    return gene_data

def getSequence(gene):
    
    EXTENSION = 6000
    loc = gene.Location
    
    #5'flanking Sequence:
    seq5Prim = ensembl.sequence.get_sequence(coord = None,
                                             genome = loc.genome,
                                             coord_name = loc.CoordName,
                                             start = loc.Start - EXTENSION,
                                             end = loc.Start,
                                             strand = 1)
        
    #gene Sequence:
    geneseq = str(gene.Seq)
                                             
    #3'flanking Sequence:
    seq3Prim = ensembl.sequence.get_sequence(coord = None,
                                             genome = loc.genome,
                                             coord_name = loc.CoordName,
                                             start = loc.End,
                                             end = loc.End + EXTENSION,
                                             strand = 1)
                                             
    return seq5Prim+geneseq+seq3Prim

##### get genes from Ensembl (pairwise or multiple)########

def get_genes_from_Ensembl_pairwise(args,sourcedata,targetdata):
    genelist = []
    
    inputtype = args.inputType
    
    firstspecies = args.firstSpecies
    secondspecies = args.secondSpecies
    if (firstspecies == None):
        print "Argument -s1 <firstspeciesname> is required"
    elif (secondspecies == None):
        print "Argument -s2 <secondspeciesname> is required"
    else:
        account = HostAccount('ensembldb.ensembl.org', 'anonymous','')
        firstgenome = Genome(firstspecies, ENSEMBL_VERSION, account)
        secondgenome = Genome(secondspecies, ENSEMBL_VERSION, account)
        
        if(inputtype == "id"):
            firstgeneid = args.firstgeneid
            if (firstgeneid == None):
                print "Argument -gid1 <firstgeneid> is required"
            secondgeneid = args.secondgeneid
            if (secondgeneid == None):
                print "Argument -gid2 <secondgeneid> is required"
            
            if(firstgeneid != None and secondgeneid != None):
                print "Retreving Genes", firstgeneid, secondgeneid

                firstgene = firstgenome.getGeneByStableId(StableId=firstgeneid)
                secondgene = secondgenome.getGeneByStableId(StableId=secondgeneid)

        if(inputtype == "name"):
            gene = args.gene
            if (gene == None):
                print "Argument -g <genename> is required"
            else:
                print "Retreving Genes", gene, "from species",  firstspecies, secondspecies

                firstgene = get_gene_from_Ensembl_by_name(gene,firstgenome)
                secondgene = get_gene_from_Ensembl_by_name(gene,secondgenome)

        sourcedata += get_cds_data(firstgene) + get_cds_data(secondgene)
        targetdata += get_gene_data(firstgene) + get_gene_data(secondgene)

    return sourcedata,targetdata

def get_genes_from_Ensembl_multiple(args,sourcedata,targetdata):
                    
    inputtype = args.inputType
    
    if(inputtype == "id"):
        geneidlistfile = args.geneIdListFile
        if (geneidlistfile == None):
            print "Argument -gidlf <geneidlistfilename> is required"
        else:
            for line in open(geneidlistfile,"r").readlines():
                parse = line.split("\n")[0].split(" ")
                if(len(parse) > 1):
                    species = parse[0]
                    geneid = parse[1]
                    print "Retreving Gene", geneid
                    account = HostAccount('ensembldb.ensembl.org', 'anonymous','')
                    genome = Genome(species, ENSEMBL_VERSION, account)
                    gene = genome.getGeneByStableId(StableId=geneid)
                    sourcedata += get_cds_data(gene)
                    targetdata += get_gene_data(gene)

    if(inputtype == "name"):
        gene = args.gene
        if (gene == None):
            print "Argument -g <genename> is required"
        specieslistfile = args.specieslistfile
        if (specieslistfile == None):
            print "Argument -slf <specieslistfilename> is required"
        if(gene != None and specieslistfile != None):
            for line in open(specieslistfile,"r").readlines():
                parse = line.split("\n")[0].split(" ")
                if(len(parse) > 0):
                    species = parse[0]
                    print "Retreving Gene", gene, "from species", species
                    account = HostAccount('ensembldb.ensembl.org', 'anonymous','')
                    genome = Genome(species, ENSEMBL_VERSION, account)
                    gene = get_gene_from_Ensembl_by_name(gene,genome)
                    sourcedata += get_cds_data(gene)
                    targetdata += get_gene_data(gene)
    return  sourcedata,targetdata

def get_gene_from_Ensembl_by_name(name,genome):
    gene = ""
    genes = genome.getGenesMatching(Symbol=name)
    for g in genes:
        gene = g
        break
    return gene


##### get data from files or Ensembl ##########

def get_data_from_files(args):
    sourcedata = []
    targetdata = []
    sourcefile = args.sourceFile
    if (sourcefile == None):
        print "Argument -sf <sourcefilename> is required"
    targetfile = args.targetFile
    if (targetfile == None):
        print "Argument -tf <targetfilename> is required"
    source2targetfile = args.source2TargetFile
    if (source2targetfile == None):
        print "Argument -s2tf <source2targetfile> is required"

    if(sourcefile != None and targetfile != None and source2targetfile != None):
        for record in SeqIO.parse(sourcefile, "fasta"):
            sourcedata.append([record.id,str(record.seq),"",[]])
        
        for record in SeqIO.parse(targetfile, "fasta"):
            targetdata.append([record.id,str(record.seq)])
        
        source2target = [line.split("\n")[0].split(" ") for line in open(source2targetfile,"r").readlines()]
        for i in range(len(sourcedata)):
            
            sourcedata[i][2] = source2target[i][1]

        sourceexonfile = args.sourceExonFile
        if(sourceexonfile != None):
            file = open(sourceexonfile, "r")
            lines = file.readlines()
            i = 0
            j = 0
            while j < len(lines):
                line = lines[j]
                if(line.startswith('>')):
                    j += 1
                    line = lines[j]
                    exonlist = []
                    while(len(line) > 2):
                        tab = line.split("\n")[0].split(" ")
                        exonlist.append([int(x) for x in tab])
                        j += 1
                        line = lines[j]
                    sourcedata[i][3] =  exonlist
                    j+=1
                    i+=1
    return  sourcedata, targetdata

def get_data_from_Ensembl(args):
    sourcedata = []
    targetdata = []
    genelist = []
    
    genenumber = args.geneNumber
    
    if(genenumber == "pairwise"):
        sourcedata,targetdata = get_genes_from_Ensembl_pairwise(args,sourcedata,targetdata)
    
    elif(genenumber == "multiple"):
        sourcedata,targetdata = get_genes_from_Ensembl_multiple(args,sourcedata,targetdata)
    
    else:
        print "Argument -gn <genenumber> is required"

    
    return sourcedata, targetdata


def get_data(args):
    sourcedata = []
    targetdata = []
    inputtype = args.inputType
    
    if (inputtype == "file"):
        sourcedata, targetdata = get_data_from_files(args)
    elif(inputtype == "id" or inputtype == "name"):
        sourcedata, targetdata = get_data_from_Ensembl(args)
    
    else:
        print "Argument -it <inputtype> is required"
    
    return sourcedata, targetdata


