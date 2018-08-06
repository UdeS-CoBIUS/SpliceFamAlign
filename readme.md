# SpliceFamAlign
Program Version 2018

compare and align all pairs of source CDS and target gene and compute transcripts orthology groups
----------------------------------------------------------------

Authors: Aïda Ouangraoua , Safa Jammali and Jean-David Aguilar
Université de Sherbrooke, Canada
Cobius Lab:  https://cobius.usherbrooke.ca/
for questions email us at Safa.Jammali@USherbrooke.ca


### Requirements:

-PyCogent
-Bio
-Ete toolkit
-Blast+ standalone
-Splign tool
-argparse


### Usage
```
usage: main.py [-h] [-c CHOICE] [-f FORCE] [-m METHOD] [-it INPUTTYPE]
               [-gn GENENUMBER] [-sf SOURCEFILE] [-tf TARGETFILE]
               [-s2tf SOURCE2TARGETFILE] [-sef SOURCEEXONFILE]
               [-s1 FIRSTSPECIES] [-s2 SECONDSPECIES] [-gid1 FIRSTGENEID]
               [-gid2 SECONDGENEID] [-gidlf GENEIDLISTFILE] [-g GENE]
               [-slf SPECIESLISTFILE] [-op OUTPUTPREFIX] [-of OUTPUTFORMAT]

```
 *-h*, --help   show this help message and exit
  
  *-c , --choice \<CHOICE>*         used where there 's no structre information, it can be blast or splign  

  *-f , --force \<FORCE>*           use SFA_G, it can be Yes or No   
  
  *-m , --method \<METHOD>*         SFA  

  *-it , --inputType \<INPUTYPE>*   file id or name   
  
  *-gn , --geneNumber \<GENENUMBER>*         Gene number: pairwise or multiple  
  
  *-sf , --sourceFile \<SOURCEFILE>*         Source file name  
  
  *-tf , --targetFile \<TARGETFILE>*         Target file name  
  
  *-s2tf , --source2TargetFile \<SOURCE2TARGETFILE>*         Source to target file name   
  
  *-sef , --sourceExonFile \<SOURCEEXONFILE>*         Source Exon file name  
  
  *-s1 , --firstSpecies \<FIRSTSPECIES>*         First species common name (required if --inputType = file or name,
                                                  and --geneNumber = pairwise")    
  
  *-s2 , -secondSpecies \<SECONDPECIES>*         Second species common name (required if --inputType = file or name,
                                                  and --geneNumber = pairwise")    
  
  *-gid1 , --firstGeneId \<FIRSTGENEID>*         First gene Ensembl Id (required if --inputType = id, and 
						--geneNumber = pairwise)    

  *-gid2 , --secondGeneId \<SECONGENEID>*         Second gene Ensembl Id (required if --inputType = id, and 
						--geneNumber = pairwise)  

  *-gidlf , --geneIdListFile \<GENEIDLISTFILE>*         Gene Id list file name (required if --inputType = id,
							 and --geneNumber = multiple)  

  *-g , -- gene \<GENE>*        Gene common name (required if --inputType = name) 
   

  *-slf , --speciesListFile \<SPECIESLISTFILE>*         Species list file name (required if --inputType = name,
                                                and --genenumber = multiple)    
  
  *-o , --outfile \<OUTFILE> *      output file   

  *-of , --outformat \<OUTFORMAT> *      output file format (list or aln)   

### Running SpliceFamAlign: example of command line


#### SFA_E when you don't have sturucture file and you use blast to compute the structure:
```
python src/main.py  -it file -sf path_to_cds_fasta_file  -tf path_to_gene_fasta_file -s2tf path_to_assiciation_cds_vs._gene_file  -op output_path -of list -m SFA -f No -c blast
```
#### SFA_E when you have you sturucture file: 
```
python src/main.py -it file -sf path_to_cds_fasta_file  -tf path_to_gene_fasta_file  -s2tf path_to_assiciation_cds_vs._gene_file -sef  path_to_structure_file -op output_path -of list -m SFA -f No 
```

 


