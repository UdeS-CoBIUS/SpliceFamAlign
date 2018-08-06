#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``write.py`` **module description**:

This module is the module that formats and writes the results in output files.

.. moduleauthor:: Aïda Ouangraoua, Safa Jammali and Jean-David Aguilar

2017-2018

"""

from alignment import *

#############################
##  INPUT FILE WRITTING #####
############################

def  write_input_files(outputprefix,sourcedata,targetdata):

    initsourcefile = open(outputprefix+"initialsource.fasta","w")
    targetfile = open(outputprefix+"target.fasta","w")
    source2targetfile = open(outputprefix+"initialsource2target.txt","w")
    initsourceexonlistfile = open(outputprefix+"initialsourceexonlist.txt","w")

    for cds in sourcedata:
        cdsid,cdsseq,cdsgeneid,exonlist = cds
	
        initsourcefile.write(">"+cdsid+"\n")
        initsourcefile.write(cdsseq+"\n\n")
        source2targetfile.write(cdsid + " " + cdsgeneid + "\n")
        initsourceexonlistfile.write(">"+cdsid+"\n")
        for block in exonlist:
            cdsstart,cdsend,genestart,geneend = block
            initsourceexonlistfile.write(str(cdsstart) + " " + str(cdsend) + " " + str(genestart) + " " + str(geneend)+"\n")
        initsourceexonlistfile.write("\n")
    initsourcefile.close()
    source2targetfile.close()
    initsourceexonlistfile.close()
    for gene in targetdata:
        geneid,geneseq = gene
        targetfile.write(">"+geneid+"\n")
        targetfile.write(geneseq+"\n\n")
    targetfile.close()

#############################
## OUTPUT FILE WRITTING #####
############################

        
def writeOutfile(outputprefix,outputformat,extendedsourcedata,targetdata,extendedorthologygroups,comparisonresults,nbinitialsource,mblocklist):
    
    sourcefile = open(outputprefix+"source.fasta","w")
    resultfile = open(outputprefix+"result.txt","w")
    source2targetfile = open(outputprefix+"source2target.txt","w")
    cdsexonendlistfile = open(outputprefix+"cdsexonendlist.txt","w")
    geneexonstartlistfile = open(outputprefix+"geneexonstartlist.txt","w")
    geneexonendlistfile = open(outputprefix+"geneexonendlist.txt","w")
    orthologygrouplistfile = open(outputprefix+"orthologygrouplist.txt","w")

    wholealignmentfile = open(outputprefix+"wholealignment.fasta","w")
    
    for cds in extendedsourcedata:
        cdsid,cdsseq,cdsgeneid ,null = cds # TO COMPUTE
        sourcefile.write(">"+cdsid+"\n")
        sourcefile.write(cdsseq+"\n\n")
        source2targetfile.write(cdsid + " " + cdsgeneid + "\n")
    sourcefile.close()
    source2targetfile.close()

    j= 0
    for gene in targetdata:
        geneid,geneseq = gene
        for i in range(nbinitialsource):
            cds = extendedsourcedata[i]
            cdsid,cdsseq,cdsgeneid ,null = cds
            cdslength = len(cdsseq)
            print "Writting results for CDS " + cdsid+ " with Gene " + geneid
            status, blocklist,splicing_sites,ttargetcds, texisting = comparisonresults[j][i]
            resultfile.write(cdsid + "\t" + geneid + "\t" + str(cdslength)  + "\t" + str("%.2f" % cover_percentage(blocklist,cdslength)) + "\t" + str(status)  + "\n")
            print_blocklist(cdsid,geneid,cdsseq, geneseq, blocklist,resultfile, outputformat)
            print_wholealignment(cdsid,geneid,cdsseq, geneseq, blocklist,wholealignmentfile, outputformat)
            
            if(cdsgeneid == geneid):
                print_exonextremitylist(cdsid, geneid, blocklist, cdsexonendlistfile,geneexonstartlistfile,geneexonendlistfile)
        j += 1
    resultfile.close()
    wholealignmentfile.close()
    cdsexonendlistfile.close()
    geneexonstartlistfile.close()
    geneexonendlistfile.close()

    if( extendedorthologygroups != []):
        for i in range(len(extendedorthologygroups)):
            orthologygrouplistfile.write(">"+str(i)+"\n")
            for k in extendedorthologygroups[i]:
                cds = extendedsourcedata[k]
                cdsid,null,null ,null = cds
                orthologygrouplistfile.write(cdsid+"\n")
            orthologygrouplistfile.write("\n")
        orthologygrouplistfile.close()

   


