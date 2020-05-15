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
(2) Fasta of nucleotide ORFs with first BLAST hit in description "{}_QMS_hits.fasta" \n\
(3) An xml of top 10 BLAST hits "{}_blast.xml" \n\
(4) Fasta of QMS transcripts which matched to the annotated transcriptome \n\
(5) Fasta of QMS transcripts which remain after failing to match the annotated transcriptome \n\
(6) New annotated transcriptome with "***" indicating which transcripts are proteomically confirmed and remaining orfs appended \n\n\
:::EXAMPLE::: \n\
QTIBAC.py -f Cline-CLP2629_combined_ORFs.fasta -q Cline-CLP2629_QMS_DB.pro -a Cline-CLP2629_transcriptome_v3.fasta -o Cline-CLP2629 -db ~/Dropbox/bin/SWISSprot_2019-09-06/SWISSprot -b blastx -p 99\n\n\
::CITE:: \n\
https://github.com/reptilerhett/Bioinformatic-Scripts\n\n')

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="REQUIRED: Nucleotide ORFs fasta made from assembly and matching the protein ORFs used for QMS identification (i.e. emboss::getorf -find 3)")
parser.add_argument("-q","--qms",
					type=argparse.FileType('r+'),
					help="\nREQUIRED: qMS protein fasta (exported from Scaffold4)")
parser.add_argument("-a","--annotated",
					type=argparse.FileType('r+'),
					help="REQUIRED: Annotated transcriptome fasta")
parser.add_argument("-db","--blast_database",
					type=str,
                    default="",
                    help="REQUIRED: Path to blast database for annotation (e.g. SWISSprot). Assumes makeblastdb has already been run.")
parser.add_argument("-b","--blast_type",
					type=str,
                    default="blastx",
                    help="Type of blast to run (e.g. blastn or blastx) against blast database for annotation. Default assumes blastx")
parser.add_argument("-c","--cdhit",
					action="store_true",
                    help="Use cd-hit-est instead of blastn to compare to your annotated transcriptome")
parser.add_argument("-o","--output",
					type=str,
                    default='QTIBAC',
                    help="Conserved part of name for output files")
parser.add_argument("-p","--matchpercent",
                    nargs='?',
                    type=float,
                    default=99,
                    help="Integer percent identity for blastn or cd-hit-est. Default is 99")
parser.add_argument("-t","--num_threads",
					type=int,
                    default=8,
                    help="Number of threads for blast. Default is 8")
parser.add_argument("-bp","--blast_path",
					type=str,
                    default="",
                    help="Directory with blast command. Default assumes it is in your path (e.g. /PATH/TO/BIN/WITH/BLAST/)")
parser.add_argument("-m","--makeblastdb",
					type=str,
                    default="makeblastdb",
                    help="Directory with makeblastdb command. Default assumes it is in your path (e.g. /PATH/TO/BIN/WITH/MAKEBLASTDB/)")
parser.add_argument("-cdp","--cdhit_path",
					type=str,
                    default="",
                    help="Directory with cd-hit-est-2d command. Default assumes it is in your path (e.g. /PATH/TO/BIN/WITH/CDHIT/)")
parser.add_argument("--version", action='version', version='QTIBAC v2')
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name
qms_name = args.qms.name
annotated_name = args.annotated.name
output = args.output
blast_database = args.blast_database
blast_path = args.blast_path
blast_type = args.blast_type
makeblastdb = args.makeblastdb
cdhit_path = args.cdhit_path
match = args.matchpercent
num_threads = args.num_threads

###############################################

print("\n::::::: Grabbing corresponding nucleotide ORFs :::::::\n")

fasta = list(SeqIO.parse(fasta_name,"fasta"))
qms = list(SeqIO.parse(qms_name,"fasta"))

names = []
for seq in qms:
	names.append(seq.id)

orfs = []
for i in names:
	for j in fasta:
		if i==j.id:
			orfs.append(j)

handle=open(output + '_QMS_DB.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(orfs)
handle.close()

###############################################

print("\n::::::: BLASTing ORFs and appending hit to description :::::::\n\n")

orfs_name = output + "_QMS_DB.fasta"
blast_xml = output + "_blast.xml"

command = blast_path + blast_type + ' -query ' + orfs_name + ' -db ' + blast_database + ' -outfmt 5 -num_threads ' + str(num_threads) + ' -max_target_seqs 10 -evalue 0.0001 -out ' + blast_xml
subprocess.call(command, shell=True)

#orfs_w_hits = []
blast_results = list(SearchIO.parse(blast_xml, 'blast-xml'))
for seq in orfs:
	for rec in blast_results:
		if seq.name == rec.id:
			hit = rec[0]
			seq.description=hit.description
			#for hit in rec:
			#	seq.description=hit.description
#				orfs_w_hits.append(seq)

handle=open(output + '_QMS_hits.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(orfs)
handle.close()

#csv_output = []
#for rec in blast_results:
#	line =[]
#	line.append(rec.id)
#	line.append(rec.seq_len)
#	for hit in rec:
#		line.append(hit.description)
#	csv_output.append(line)

#with open(output + '_hits.csv', 'w') as csv_file:
#	csv_writer = csv.writer(csv_file, delimiter = ',')
#	for row in csv_output:
#		csv_writer.writerow(row)

###############################################

print("\n::::::: Checking against annotated transcriptome :::::::\n")

annotated = list(SeqIO.parse(annotated_name,"fasta"))

if args.cdhit:
	match = match/100
	command = cdhit_path + "cd-hit-est-2d -i " + annotated_name + " -i2 " + orfs_name + " -o tmp.fasta -d 0 -c " + str(match)
	subprocess.call(command,shell=True)
	
	clustered = open('tmp.fasta.clstr',"r")
	clustered = clustered.read()
	clustered = clustered.split('>Cluster ')[1:]
	
	clustered_array = []
	for line in clustered:
		clustered_array.append(line.split('\n')[:-1])
	
	confirmed_transcripts = []
	confirmed_orfs = []
	confirmed_orfs_names = []
	for cluster in clustered_array:
		if len(cluster) > 2:
			for i in range(1,len(cluster)):
				name = cluster[i].split(', >')[1].split('...')[0]
				for seq in annotated:
					if name == seq.name :
						seq.id=seq.name + "***"
						seq.description=seq.name + "***"
						seq.name=seq.id
						confirmed_transcripts.append(seq.name)
				for seq in orfs:
					if name == seq.name :
						confirmed_orfs.append(seq)
						confirmed_orfs_names.append(seq.name)
	
	remaining_orfs = [x for x in orfs if x.name not in confirmed_orfs_names]
	
	handle=open(output + '_QMS_matched_orfs.fasta', "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(confirmed_orfs)
	handle.close()
	
	handle=open(output + '_QMS_remaining_orfs.fasta', "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(remaining_orfs)
	handle.close()
	
	command = cdhit_path + "cd-hit-est -i " + output + "_QMS_remaining_orfs.fasta -o " + output + "_QMS_remaining_orfs_clustered.fasta -d 0 -c " + str(match)
	subprocess.call(command,shell=True)
	
	clustered_remaining_orfs = list(SeqIO.parse(output + "_QMS_remaining_orfs_clustered.fasta","fasta"))
	
	new_annotated = annotated.copy()
	for seq in clustered_remaining_orfs:
		seq.id="PUTATIVE_" + seq.name + "***"
		new_annotated.append(seq)
	
	handle=open(output + '_new_transcriptome.fasta', "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(new_annotated)
	handle.close()
	
	print("\n\n::::::: RESULTS :::::::\n")
	print("Used cd-hit to determine that:")
	print(str(len(confirmed_orfs)) + " of the " + str(len(orfs)) + " QMS-identified ORFs match previously annotated transcripts")
	print("These matches are equal to " + str(len(confirmed_transcripts)) + " proteomically-confirmed transcripts: \n")
	
	print(*confirmed_transcripts, sep = "\n") 
	
	print("\n\n" + str(len(remaining_orfs)) + " of the " + str(len(orfs)) + " QMS-identified ORFs did not match previously annotated transcripts")
	print("After clustering, this resulted in " + str(len(clustered_remaining_orfs)) + " QMS-identified ORFs which have been added to the end of your annotated transcriptome")
	
	subprocess.call('rm tmp.fasta*',shell=True)


else:
	command = makeblastdb + " -dbtype nucl -in " + annotated_name + " -out tmp.blastdb"
	subprocess.call(command,shell=True)

	command =  blast_path + "blastn -query " + orfs_name + " -out match.blast.out -db tmp.blastdb -perc_identity " + str(match) + " -outfmt 5 -num_threads " + str(num_threads) + " -evalue 0.0001"
	subprocess.call(command,shell=True)

	blast_results = list(SearchIO.parse("match.blast.out", 'blast-xml'))

	confirmed_orfs = []
	remaining_orfs = []
	confirmed_orfs_names = []
	remaining_orfs_names = []
	for seq in orfs:
		for rec in blast_results:
			if seq.name == rec.id:
				if len(rec.hits) == 0:
					remaining_orfs.append(seq)
					remaining_orfs_names.append(seq.name)
				else:
					#seq.description = rec.hits[0].id
					confirmed_orfs.append(seq)
					confirmed_orfs_names.append(seq.name)

	handle=open(output + '_QMS_matched_orfs.fasta', "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(confirmed_orfs)
	handle.close()
	
	handle=open(output + '_QMS_remaining_orfs.fasta', "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(remaining_orfs)
	handle.close()

	command = cdhit_path + "cd-hit-est -i " + output + "_QMS_remaining_orfs.fasta -o " + output + "_QMS_remaining_orfs_clustered.fasta -d 0 -c " + str(match/100)
	subprocess.call(command,shell=True)

	clustered_remaining_orfs = list(SeqIO.parse(output + "_QMS_remaining_orfs_clustered.fasta","fasta"))

	confirmed_transcripts = []
	for seq in annotated:
		for rec in blast_results:
			for hit in rec:
				if seq.name==hit.id:
					seq.id=seq.name + "***"
					seq.description=seq.name + "***"
					seq.name=seq.id
					if seq.name not in confirmed_transcripts:
						confirmed_transcripts.append(seq.name)

	new_annotated = annotated.copy()
	for seq in clustered_remaining_orfs:
		seq.id="PUTATIVE_" + seq.name + "***"
		new_annotated.append(seq)

	handle=open(output + '_new_transcriptome.fasta', "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(new_annotated)
	handle.close()

	csv_output = []
	for rec in blast_results:
		line =[]
		line.append(rec.id)
		line.append(rec.seq_len)
		for hit in rec:
			line.append(hit.id)
		csv_output.append(line)

	with open(output + '_QMS_matches.csv', 'w') as csv_file:
		csv_writer = csv.writer(csv_file, delimiter = ',')
		for row in csv_output:
			csv_writer.writerow(row)

	print("\n\n::::::: RESULTS :::::::\n")
	print("Used blast to determine that:")
	print(str(len(confirmed_orfs)) + " of the " + str(len(orfs)) + " QMS-identified ORFs matched previously annotated transcripts")
	print("These matches are equal to " + str(len(confirmed_transcripts)) + " proteomically-confirmed transcripts: \n")

	print(*confirmed_transcripts, sep = "\n") 

	print("\n\n" + str(len(remaining_orfs)) + " of the " + str(len(orfs)) + " QMS-identified ORFs did not match previously annotated transcripts")
	print("After clustering, this resulted in " + str(len(clustered_remaining_orfs)) + " QMS-identified ORFs which have been added to the end of your annotated transcriptome")

	subprocess.call('rm match.blast.out',shell=True)
	subprocess.call('rm tmp.blastdb*',shell=True)
	subprocess.call('rm ' + blast_xml, shell=True)

###############################################

print("\n::::::: FINISHED :::::::\n")