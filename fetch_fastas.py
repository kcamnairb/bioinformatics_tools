import argparse
import subprocess,re, sys
from os import listdir
from os.path import isfile, join
from Bio import AlignIO, SeqIO
from io import StringIO
parser = argparse.ArgumentParser(description="Retrieves individual fastas from multifasta files. Each fasta ID that you want to fetch should be on its own line.")
parser.add_argument("fastas_to_fetch")
parser.add_argument("fasta_files", nargs='+')
args = parser.parse_args()

fasta_dict = {}
for file in args.fasta_files:
    fasta_dict.update(SeqIO.to_dict(SeqIO.parse(open(file),"fasta")))
with open(args.fastas_to_fetch) if args.fastas_to_fetch is not "-" else sys.stdin as f:
	target_fastas = [fasta.strip().split(" ")[0] for fasta in f]
num_fastas_fetched = 0
missing = []
with open("fetched_fastas.fasta","w") as outfile:
	for fasta in target_fastas:
		if fasta in fasta_dict:
			SeqIO.write(fasta_dict[fasta],outfile,"fasta")
			num_fastas_fetched += 1
		else:
			missing.append(fasta)
print(str(num_fastas_fetched) + ' fastas fetched')
if missing:
	print('Could not find ' + str(', '.join(missing)))