#!/usr/bin/env python
import argparse
import os.path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Calculates N50 from one or more fasta files, and makes a histogram")
parser.add_argument("fasta_files", nargs='+')
args = parser.parse_args()
fasta_files = args.fasta_files

def calculate_n50(fasta):
    seq_lengths = sorted([len(seq) for seq in SeqIO.parse(open(fasta), "fasta")], reverse=True)
    seq_lengths_above_1000 = [length for length in seq_lengths if length >= 1000]
    genome_size = sum(seq_lengths_above_1000) 
    half_genome = genome_size / 2
    test_sum = genome_size
    for contig_num, length in enumerate(seq_lengths_above_1000):
        test_sum -= length
        if test_sum < half_genome:
            N50_contig = contig_num - 1
            N50 = seq_lengths_above_1000[N50_contig]
            break
    return (([os.path.basename(fasta), str(sum(seq_lengths)), str(sum(seq_lengths_above_1000)), str(N50), \
        str(N50_contig), str(len(seq_lengths)), str(len(seq_lengths_above_1000)), str(seq_lengths[0])]), seq_lengths_above_1000)

with open('genome_stats.csv','w') as out_f:
    out_f.write(','+ ','.join(['total_length', 'total_length_of_contigs_above_1000',  'N50', \
        'N50_contig_num', 'number_of_contigs', 'number_of_contigs_above_1000', 'max_contig_size']) +'\n')    
    stats = [calculate_n50(fasta) for fasta in fasta_files]
    for fasta in stats:
        out_f.write(','.join(fasta[0])+'\n')
titles = [assembly[0][0] for assembly in stats]
data = [assembly[1] for assembly in stats]
if len(fasta_files) > 1:
    fig, ax = plt.subplots(len(data), 1, sharex=True, sharey=True)
    ax = ax.ravel()
    for i, ax in enumerate(ax):
        ax.hist(data[i], bins=25)
        ax.set_title(titles[i])
else:
    plt.hist(data[0], bins=25)
plt.ylabel('# of contigs')
plt.xlabel('length of contig')
plt.tight_layout()
plt.savefig('histogram.png')
