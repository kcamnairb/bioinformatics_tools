import os
from flask import Flask, request, redirect, url_for, flash
from werkzeug import secure_filename
import subprocess
from os.path import isfile, basename, splitext
from Bio.SeqUtils import GC
from Bio import SeqIO
from flask import send_from_directory
# Takes any combination of nucletide or amino acid fasta files, makes a database of the subject file, and
# determines the most appropriate type of BLAST program to run and executes it with the given query and subject files. 
# Requires BLAST+, Biopython and Flask to be installed
UPLOAD_FOLDER = '.'
ALLOWED_EXTENSIONS = set(['txt','fasta','fa', 'fsa', 'fna', 'fas'])
os.chdir('/home/brian/flask')
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
	return '.' in filename and \
		   filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
	if request.method == 'POST':
		query_file = request.files['query']
		if query_file and allowed_file(query_file.filename):
			query_filename = secure_filename(query_file.filename)
			query_file.save(os.path.join(app.config['UPLOAD_FOLDER'], query_filename))
		subject_file = request.files['subject']
		if subject_file and allowed_file(subject_file.filename):
			subject_filename = secure_filename(subject_file.filename)
			subject_file.save(os.path.join(app.config['UPLOAD_FOLDER'], subject_filename))
		if request.form['submit'] == 'BLAST' and query_file and subject_file:
			blast(query_filename, subject_filename)
			return redirect(url_for('uploaded_file', 
			filename=splitext(basename(query_filename))[0] + '_vs_' + splitext(basename(subject_filename))[0]))

	return '''
	<!doctype html>
	<title>BLAST two multifasta files</title>
	<h2>Upload Query File</h2>
	<form action="" method=post enctype=multipart/form-data>
	  <p><input type=file name=query>
	<h2>Upload Subject File</h2>
	<form action="" method=post enctype=multipart/form-data>
	  <p><input type=file name=subject>
	<br/>
	<br/>
	<input type=submit name=submit value=BLAST>
	</form>
	'''
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

	subprocess.call('makeblastdb -parse_seqids -dbtype ' + dbtype + ' -in ' + subject_file + ' -out ' + subject_file, shell=True) 
	subprocess.call (blast_program + ' -query ' + query_file + ' -db ' + subject_file + ' -out ' + splitext(basename(query_file))[0] + '_vs_' + splitext(basename(subject_file))[0] +
		' -num_threads 16 -evalue 1e-10 -max_target_seqs 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qcovs" ',shell=True)
	outfile = splitext(basename(query_file))[0] + '_vs_' + splitext(basename(subject_file))[0]
	temp = open(outfile).read()
	results_w_header = open(outfile,'w')
	header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle', 'qcovs']	
	results_w_header.write('\t'.join(header)+'\n')
	for line in temp:
		results_w_header.write(line)
	results_w_header.close
	

@app.route('/uploads/<filename>')
def uploaded_file(filename):
	return send_from_directory(app.config['UPLOAD_FOLDER'],
							   filename)

if __name__ == '__main__':
	app.run(host='0.0.0.0', debug = False)