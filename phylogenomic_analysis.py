#!/usr/bin/env python
import subprocess 
import re
import sys
import argparse
import os
from os.path import isfile, join, isdir
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO, SeqIO
from io import StringIO
from collections import namedtuple

parser = argparse.ArgumentParser(description='Aligns orthologous protein groups and creates phylogenetic tree. Requires Muscle and Biopython')
parser.add_argument('-orth', '--ortholog_file', help='Proteinortho output file or tab separated list of orthologs', required=True)
parser.add_argument('-outg', '--outgroup', help='Name of outgroup, should match one of species from header of ortholog file', required=True)
parser.add_argument('-fasta', '--fasta_directory', help='Directory containing fasta files with sequences of all proteins in ortholog file', required=True)
args = parser.parse_args()
fasta_dir = args.fasta_directory
fasta_files = [join(fasta_dir,f) for f in os.listdir(fasta_dir) if isfile(join(fasta_dir,f))] #create list of fasta files from fasta_file directory
protein_list_raw = [protein_group.strip().split('\t') for protein_group in open(args.ortholog_file,'rU') if len(protein_group.split(',')) < 2 and '*' not in protein_group]
Protein = namedtuple('Protein', ['name', 'species'])
temp = []
if '# Species' in protein_list_raw[0][0]:
    protein_list_raw = [protein_group[3:] for protein_group in protein_list_raw ]
else:
    print('"# Species" not found, assuming not proteinortho format, using all columns')
header = protein_list_raw.pop(0)
num_species = len(header)
protein_list = [[Protein(*protein) for protein in zip(protein_group, header)] for protein_group in protein_list_raw]
#paralogs are seperated by "," so any ortholog cluster containing paralogs is not used
fasta_objects = []
for multifasta in fasta_files:
    fasta_objects.extend(SeqIO.parse(multifasta,"fasta"))
protein_fasta_dict = SeqIO.to_dict(fasta_objects) 
protein_file_list = []
for idx, _ in enumerate(protein_list):
    protein_file_list.append("clusters/protein_group"+str(idx)+".fasta") #creates a list of fasta files to write to in the next step
if not isdir('clusters'):
    os.makedirs('clusters')
for idx, cluster in enumerate(protein_list):
    with open(protein_file_list[idx],"w") as outfile:
        for protein in cluster:
            if protein.name in protein_fasta_dict:
                protein_record = protein_fasta_dict[protein.name]
                protein_record.id = protein.species+'|'+protein.name #prepends the species name to protein in case the protein file doesn't have consistent prefixes
                SeqIO.write(protein_record, outfile, "fasta")
            else:
                print('failed to find ' + protein)
                sys.exit(0)
                
alignment_file_list = []
alignment_file_list = [protein_file + '.align' for protein_file in protein_file_list] #create list of alignment files

for idx, alignfile in enumerate(alignment_file_list):
    with open(alignfile,"w") as outfile:
        cline = MuscleCommandline(input=protein_file_list[idx])
        stdout, stderr = cline() #performs muscle alignment
        align = AlignIO.read(StringIO(stdout), "fasta") #creates an alignment record
        align.sort() #sorts alignment alphabetically so that each species will be in the same position for Gblocks concatenation
        AlignIO.write(align,outfile,"fasta")
subprocess.call('ls clusters/protein_group*.fasta | parallel --jobs 31 muscle -in {} -out {}.align', shell=True)
for file in alignment_file_list:
    with open(file, 'r+') as f:
        align = AlignIO.read(f, "fasta")
        align.sort()
        f.seek(0)
        AlignIO.write(align,f,"fasta")
with open("alignment_paths", "w") as path_file:
    for file in alignment_file_list:
        path_file.write(file+"\n") #creates a file containing the paths of the alignment files for Gblocks to use in next step
subprocess.call('sed -ri "s/\|.*//" clusters/*', shell=True)
subprocess.call(["Gblocks", "alignment_paths", "-t=p", "-a=y"])
subprocess.call(["raxmlHPC-PTHREADS-SSE3","-s","clusters/alignment_paths-gb.seq","-N","autoMRE",\
    "-n","result","-k","-f","a","-p","12345","-x","12345","-m","PROTCATRTREV","-o",args.outgroup,"-T","31"])