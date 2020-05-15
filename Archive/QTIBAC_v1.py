#!/usr/bin/env python

import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='QTIBAC: QMS-based Transcript Identification, Blast, Annotation, and Confirmation \n\n\
This script is designed to deal with proteomic annotation and confirmation. \n\n\
(1) QTIBAC will take a fasta of nucleotide ORFs identified from your assembly (i.e. emboss::getorf -find 3) and qMS protein sequences (exported from Scaffold4) to identify the original nucleotide sequences. \n\
(2) QTIBAC will blast these sequences against a specified database and append the name of the first hit. \n\
(3) Finally, QTIBAC will compare the QMS transcripts to your annotated transcriptome to determine which transcripts are proteomically confirmed and what is potentially new (no matches). \n\n\
:::OUTPUT::: \n\
(1) Fasta of nucleotide ORFs matching QMS proteins "{}_QMS_DB.fasta" \n\
(2) Fasta of nucleotide ORFs with database BLAST hits "{}_QMS_hits.fasta" \n\
(3) csv of BLAST hits to annotated transcriptome "{}_hits.csv" \n\
(4) Fasta of QMS transcripts without matches in the annotated transcriptome \n\
(5) Fasta of QMS transcripts with matches in the annotated transcriptome \n\
(6) Annotated transcriptome with "***" indicating which transcripts are proteomically confirmed \n\n\
:::EXAMPLE::: \n\
QTIBAC.py -f Cline-CLP2629_combined_ORFs.fasta -q Cline-CLP2629_QMS_DB.pro -a Cline-CLP2629_transcriptome_v3.fasta -o Cline-CLP2629 -db ~/Dropbox/bin/SWISSprot_2019-09-06/SWISSprot -b blastx \n\n')
parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="REQUIRED: Nucleotide ORFs fasta made from assembly and matching the protein ORFs used for QMS identification (i.e. emboss::getorf -find 3)")
parser.add_argument("-q","--qms",
					type=argparse.FileType('r+'),
					help="REQUIRED: qMS protein fasta (exported from Scaffold4)")
parser.add_argument("-a","--annotated",
					type=argparse.FileType('r+'),
					help="REQUIRED: Annotated transcriptome fasta")
parser.add_argument("-o","--output",
					type=str,
                    default='',
                    help="Conserved part of name for output files")
parser.add_argument("-db","--blast_database",
					type=str,
                    default="",
                    help="Path to blast database for annotation (e.g. SWISSprot)")
parser.add_argument("-b","--blast_type",
					type=str,
                    default="blastn",
                    help="Type of blast to run (e.g. blastn or blastx) against blast database for annotation. Default assumes blastn")
parser.add_argument("-bp","--blast_path",
					type=str,
                    default="",
                    help="Directory with blast command. Default assumes it is in your path (e.g. /PATH/TO/BIN/WITH/BLAST/)")
parser.add_argument("-m","--makeblastdb",
					type=str,
                    default="makeblastdb",
                    help="Directory with makeblastdb command. Default assumes it is in your path")
parser.add_argument("-p","--matchpercent",
                    nargs='?',
                    type=float,
                    default=100,
                    help="Percent identity for blastn. Default is 100")
parser.add_argument("-nt","--num_threads",
					type=int,
                    default=8,
                    help="Number of threads for blast. Default is 8")
args=parser.parse_args()

###############################################

fasta = args.fasta
qms = args.qms
annotated = args.annotated
annotated_name = args.annotated.name
output = args.output
blast_database = args.blast_database
blast_path = args.blast_path
blast_type = args.blast_type
makeblastdb = args.makeblastdb
match = args.matchpercent
num_threads = args.num_threads

###############################################

print("\n::::::: Grabbing corresponding nucleotide ORFs :::::::\n")

nucleotides = list(SeqIO.parse(fasta,"fasta"))
proteins = list(SeqIO.parse(qms,"fasta"))

names = []
for seq in proteins:
	names.append(seq.id)

orfs = []
for i in names:
	for j in nucleotides:
		if i==j.id:
			orfs.append(j)

#for seq in orfs:
#	seq.description=""


handle=open(output + '_QMS_DB.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(orfs)
handle.close()

###############################################

print("\n::::::: BLASTing ORFs and appending hit to description :::::::\n\n v Ignore this warning v \n")

orfs_name = output + "_QMS_DB.fasta"
blast_xml = output + "_blast.xml"

command = blast_path + blast_type + ' -query ' + orfs_name + ' -db ' + blast_database + ' -outfmt 5 -num_threads ' + str(num_threads) + ' -max_target_seqs 1 -evalue 0.0001 -out ' + blast_xml
subprocess.call(command, shell=True)

orfs_w_hits = []
blast_results = list(SearchIO.parse(blast_xml, 'blast-xml'))
for seq in orfs:
	for rec in blast_results:
		if seq.name == rec.id:
			for hit in rec:
				seq.description=hit.description
				orfs_w_hits.append(seq)

handle=open(output + '_QMS_hits.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(orfs_w_hits)
handle.close()

###############################################

print("\n::::::: Checking against annotated transcriptome :::::::\n")

annotated = list(SeqIO.parse(annotated,"fasta"))

command = makeblastdb + " -dbtype nucl -in " + annotated_name + " -out tmp.blastdb"
subprocess.call(command,shell=True)

command =  blast_path + "blastn -query " + orfs_name + " -out match.blast.out -db tmp.blastdb -perc_identity " + str(match) + " -outfmt 5 -num_threads " + str(num_threads) + " -evalue 0.0001"
subprocess.call(command,shell=True)

blast_results = list(SearchIO.parse("match.blast.out", 'blast-xml'))

qms_no_match = []
qms_match = []
for seq in orfs_w_hits:
	for rec in blast_results:
		if seq.name == rec.id:
			if len(rec.hits) == 0:
				qms_no_match.append(seq)
			else:
				#seq.description = rec.hits[0].id
				qms_match.append(seq)

for seq in annotated:
	for rec in blast_results:
		for hit in rec:
			if seq.name==hit.id:
				seq.id=seq.name + "***"
				break

count = []
for seq in annotated:
	if "***" in seq.id:
		count.append(1)

csv_output = []
for rec in blast_results:
	line =[]
	line.append(rec.id)
	line.append(rec.seq_len)
	for hit in rec:
		line.append(hit.id)
	csv_output.append(line)

with open(output + '_hits.csv', 'w') as csv_file:
	csv_writer = csv.writer(csv_file, delimiter = ',')
	for row in csv_output:
		csv_writer.writerow(row)

print("\n\n" + str(len(qms_match)) + " of the " + str(len(orfs_w_hits)) + " QMS-identified ORFs match previously annotated transcripts")
print("These matches are equal to " + str(sum(count)) + " proteomically-confirmed transcripts")

print("\n\n" + str(len(qms_no_match)) + " of the " + str(len(orfs_w_hits)) + " QMS-identified ORFs with no matches.")
print("You may want to check these and add them to your transcriptome?")

###############################################

handle=open(output + '_QMS_transcriptome_no_match.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(qms_no_match)
handle.close()

handle=open(output + '_QMS_transcriptome_match.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(qms_match)
handle.close()

handle=open(output + '_transcriptome_QMS_confirmed.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(annotated)
handle.close()

subprocess.call("perl -pi -e 's/ .*//g' *_transcriptome_QMS_confirmed.fasta", shell=True)
subprocess.call('rm match.blast.out',shell=True)
subprocess.call('rm tmp.blastdb*',shell=True)
subprocess.call('rm ' + blast_xml, shell=True)

###############################################

print("\n::::::: FINISHED :::::::\n")