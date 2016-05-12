import os, argparse
import subprocess
from os.path import isfile, basename, splitext
from Bio.SeqUtils import GC
from Bio import SeqIO
import multiprocessing


def blast(query, subject):

	query_file = query
	subject_file = subject
	subject_record_list = [record for record in list(SeqIO.parse(subject_file, "fasta"))]
	query_record_list = [record for record in list(SeqIO.parse(query_file, "fasta"))]
	subject_length = sum([len(rec) for rec in subject_record_list])
	query_length = sum([len(rec) for rec in query_record_list])
	subject_gc = sum([GC(rec.seq)*len(rec) for rec in subject_record_list])/subject_length
	query_gc = sum([GC(rec.seq)*len(rec) for rec in query_record_list])/query_length

	dbtype, blast_program = '', ''
	if subject_gc > 25 and query_gc > 25:
		dbtype, blast_program = 'nucl', 'blastn'
	elif subject_gc < 25 and query_gc < 25:
		dbtype, blast_program = 'prot', 'blastp'
	elif subject_gc > 25 and query_gc < 25:
		dbtype, blast_program = 'nucl', 'tblastn'
	elif subject_gc < 25 and query_gc > 25:
		dbtype, blast_program = 'prot', 'blastx'
    num_threads = str(multiprocessing.cpu_count() - 1)
	subprocess.call('makeblastdb -parse_seqids -dbtype ' + dbtype + ' -in ' + subject_file + ' -out ' + subject_file, shell=True) 
	subprocess.call (blast_program + ' -query ' + query_file + ' -db ' + subject_file + ' -out ' + splitext(basename(query_file))[0] + '_vs_' + splitext(basename(subject_file))[0] +
		' -num_threads '+num_threads+ ' -evalue 1e-10 -max_target_seqs 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qcovs" ',shell=True)
	outfile = splitext(basename(query_file))[0] + '_vs_' + splitext(basename(subject_file))[0]
	temp = open(outfile).read()
	results_w_header = open(outfile,'w')
	header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle', 'qcovs']
	results_w_header.write('\t'.join(header)+'\n')
	for line in temp:
		results_w_header.write(line)
	results_w_header.close
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='''Takes any combination of nucletide or amino acid fasta files, makes a database of the subject file, and
    determines the most appropriate type of BLAST program to run and executes it with the given query and subject files. Requires BLAST+ and Biopython to be installed.''')
	parser.add_argument('-query', '--query_file', required=True)
	parser.add_argument('-subject', '--subject_file', required=True)
	args = parser.parse_args()
	blast(args.query_file, args.subject_file)