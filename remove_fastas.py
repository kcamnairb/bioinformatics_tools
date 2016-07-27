#!/usr/bin/env python3
import argparse
import os
from Bio import SeqIO
parser = argparse.ArgumentParser(description="Remove individual fastas from multifasta file. Each fasta ID that you want to remove should be on its own line.")
parser.add_argument("id_list")
parser.add_argument("fasta_file")
args = parser.parse_args()
id_list =  os.path.splitext(os.path.basename(args.id_list))[0]
outfasta, ext = os.path.splitext(os.path.basename(args.fasta_file))
outfasta = outfasta + '_' + id_list + '_removed.' + ext

fasta_dict = SeqIO.to_dict(SeqIO.parse(open(args.fasta_file),"fasta"))
target_fastas = [fasta.strip().split(" ")[0] for fasta in open(args.id_list)]
removed = []
with open(outfasta,"w") as outfile:
    for fasta in fasta_dict:
        if fasta not in target_fastas:
            SeqIO.write(fasta_dict[fasta],outfile,"fasta")
        else:
            removed.append(fasta)
    print(str(len(target_fastas)) + ' asked to remove, ' + str(len(removed)) + ' removed')