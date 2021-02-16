#!/home3/sthio/anaconda/envs/py356/bin/python3.5
#Author: Shinta Thio
#Last updated: Dec 3, 2019
#Purpose: scan for upstream and downstream start codon from the strast codon predicted by transdecoder

import os, sys, pyfaidx, re
import gffutils as gff
from gffutils import Feature

translation_block_gff3 = sys.argv[1]
qry		 			   = sys.argv[2]

qry_db_name = qry + '.db'
fasta = pyfaidx.Fasta('/home3/sthio/PROJECTS/C_briggsae101/gene_model/c_briggsae.PRJNA10731.WS254.genomic.fa')


if os.path.isfile(qry_db_name):
    print('Query database found. Importing...', file=sys.stderr)
    qry_db = gff.FeatureDB(qry_db_name)
else:
    print('Building query database...', file=sys.stderr)
    qry_db = gff.create_db(qry, dbfn=qry_db_name, 
        merge_strategy='create_unique', 
        disable_infer_genes=True, 
        disable_infer_transcripts=True,
        verbose=True,
        keep_order=True)

def correct_startcodon_plus (q_transcript, first_cds_seq, first_cds_s, first_cds_e):
	first_cds_seq_candidate = ''
	corrected_first_cds_seq = ''
	corrected_first_seq_s = 0
	atg_location = ''
	length = len(first_cds_seq)
	plus_dicts=[plus0_dict, plus1_dict, plus2_dict]
	for plus_dict in plus_dicts: 
		for coordinate, block_seq in plus_dict.items():
			if first_cds_seq in str(block_seq):
				if len(block_seq)%3 == 0:
					for i in range(0, len(block_seq), 3):
						seq_candidate = block_seq[i:len(block_seq)]
						if str(seq_candidate).startswith('ATG'):
							first_cds_seq_candidate = ''.join(str(seq_candidate).partition(first_cds_seq)[0:2])
							if first_cds_seq_candidate.endswith(first_cds_seq):    								#ATG upstream
								corrected_first_cds_seq = first_cds_seq_candidate
								corrected_first_seq_s = first_cds_s - (len(corrected_first_cds_seq)-len(first_cds_seq))
								atg_location = 'upstream'
							else: 																						#ATG downstream
								if length >= 9:
									first_cds_seq_candidate_partition = first_cds_seq_candidate.split((first_cds_seq[-9:]))
									corrected_first_cds_seq = first_cds_seq_candidate_partition[0]+(first_cds_seq[-9:])
									corrected_first_seq_s = first_cds_s + (len(first_cds_seq)-len(corrected_first_cds_seq))
									atg_location = 'downstream'
								else:
									first_cds_seq_candidate_partition = first_cds_seq_candidate.split((first_cds_seq[-(length):]))
									corrected_first_cds_seq = first_cds_seq_candidate_partition[0]+(first_cds_seq[-(length):])
									corrected_first_seq_s = first_cds_s + (len(first_cds_seq)-len(corrected_first_cds_seq))
									atg_location = 'downstream'
							break
	return corrected_first_cds_seq, corrected_first_seq_s, atg_location

def correct_startcodon_minus (q_transcript,first_cds_seq, first_cds_s, first_cds_e):
	first_cds_seq_candidate = ''
	corrected_first_cds_seq = ''
	corrected_first_seq_e = 0
	atg_location = ''
	length = len(first_cds_seq)
	minus_dicts=[minus0_dict, minus1_dict, minus2_dict]
	for minus_dict in minus_dicts:
		for coordinate, block_seq in minus_dict.items():
			if first_cds_seq in str(block_seq):
				if len(block_seq)%3 == 0:
					for i in range(0, len(block_seq), 3):
						seq_candidate = block_seq[i:len(block_seq)]
						if str(seq_candidate).startswith('ATG'):
							first_cds_seq_candidate = ''.join(str(seq_candidate).partition(first_cds_seq)[0:2])
							if first_cds_seq_candidate.endswith(first_cds_seq): 									#ATG upstream
								corrected_first_cds_seq = first_cds_seq_candidate
								corrected_first_seq_e = first_cds_e + (len(corrected_first_cds_seq)-len(first_cds_seq))
								atg_location = 'upstream'
							else: 																						#ATG downstream
								if length >= 9:
									first_cds_seq_candidate_partition = first_cds_seq_candidate.split((first_cds_seq[-9:]))
									corrected_first_cds_seq = first_cds_seq_candidate_partition[0]+(first_cds_seq[-9:])
									corrected_first_seq_e = first_cds_e - (len(first_cds_seq)-len(corrected_first_cds_seq))
									atg_location = 'downstream'
								else:
									first_cds_seq_candidate_partition = first_cds_seq_candidate.split((first_cds_seq[-(length):]))
									corrected_first_cds_seq = first_cds_seq_candidate_partition[0]+(first_cds_seq[-(length):])
									corrected_first_seq_e = first_cds_e - (len(first_cds_seq)-len(corrected_first_cds_seq))
									atg_location = 'downstream'
							break
	return corrected_first_cds_seq, corrected_first_seq_e, atg_location

def output_features(q_transcript):
	genes = qry_db.parents(q_transcript, featuretype='gene', order_by='start')
	for gene in genes:
		print(gene)
	print(q_transcript)
	for feature in qry_db.children(q_transcript, featuretype=['five_prime_UTR','exon','CDS','intron','three_prime_UTR'], order_by='start'): 
		print(feature)

def output_features_plus_upstream(q_transcript, corrected_first_seq_s):
	genes = qry_db.parents(q_transcript, featuretype='gene', order_by='start')
	for gene in genes:
		gene.start = corrected_first_seq_s
		print(gene)
	q_transcript.start=corrected_first_seq_s
	print(q_transcript)
	for feature in qry_db.children(q_transcript, featuretype=['five_prime_UTR','exon','CDS','intron','three_prime_UTR'], order_by='start'): 
		if feature[2] == 'CDS':
			if first_cds_s == feature.start:
				feature.start = corrected_first_seq_s
				print(feature)
			else:
				print(feature)
		if feature[2] == 'exon':
			if first_cds_s == feature.start:
				feature.start = corrected_first_seq_s
				print(feature)
			else:
				print(feature)
		if feature[2]!='CDS' and feature[2]!='exon':
			print(feature)

def output_features_plus_downstream(q_transcript, corrected_first_seq_s):
	genes = qry_db.parents(q_transcript, featuretype='gene', order_by='start')
	for gene in genes:
		print(gene)
	print(q_transcript)
	for feature in qry_db.children(q_transcript, featuretype=['five_prime_UTR','exon','CDS','intron','three_prime_UTR'], order_by='start'): 
		if feature.featuretype == 'CDS':
			if first_cds_s == feature.start:
				feature.start = corrected_first_seq_s
				if not 'XLOC' in q_transcript_id:
					utr_attr = 'Parent=Transcript:'+tracking_q_transcript_id
					line = q_chrom, 'transdecoder', 'five_prime_UTR', first_cds_s, (corrected_first_seq_s-1), '.', q_strand, '.', utr_attr
					print(*line, sep='\t')
				else:
					utr_attr = 'Parent='+transcript_id
					line = q_chrom, 'transdecoder', 'five_prime_UTR', first_cds_s, (corrected_first_seq_s-1),  '.', q_strand, '.', utr_attr
					print(*line, sep='\t')
				print(feature)
			else:
				print(feature)
		if feature.featuretype!='CDS':
			print(feature)

def output_features_minus_upstream(q_transcript, corrected_first_seq_e):#, f):
	genes = qry_db.parents(q_transcript, featuretype='gene', order_by='start')
	for gene in genes:
		gene.end = corrected_first_seq_e
		print(gene)
	q_transcript.end = corrected_first_seq_e
	print(q_transcript)
	for feature in qry_db.children(q_transcript, featuretype=['five_prime_UTR','exon','CDS','intron','three_prime_UTR'], order_by='start'): 
		#print(feature)
		if feature[2]=='CDS':
			if first_cds_e == feature.end:
				feature.end = corrected_first_seq_e
				print(feature)
			else:
				print(feature)
		if feature[2]=='exon':
			if first_cds_s == feature.start:
				feature.end = corrected_first_seq_e
				print(feature)
			else:
				print(feature)
		if feature[2]!='CDS' and feature[2]!='exon':
			print(feature)

def output_features_minus_downstream(q_transcript, corrected_first_seq_e):#, f):
	genes = qry_db.parents(q_transcript, featuretype='gene', order_by='start')
	for gene in genes:
		print(gene)
	print(q_transcript)
	for feature in qry_db.children(q_transcript, featuretype=['five_prime_UTR','exon','CDS','intron','three_prime_UTR'], order_by='start'): 
		if feature[2]=='CDS':
			if first_cds_e == feature.end:
				feature.end = corrected_first_seq_e
				print(feature)
				if not 'XLOC' in q_transcript_id:
					utr_attr = 'Parent=Transcript:'+tracking_q_transcript_id
					line = q_chrom, 'transdecoder', 'five_prime_UTR', (corrected_first_seq_e+1), q_transcript.end, '.', q_strand, '.', utr_attr
					print(*line, sep='\t')
				else:
					utr_attr = 'Parent='+transcript_id
					line = q_chrom, 'transdecoder', 'five_prime_UTR', (corrected_first_seq_e+1), q_transcript.end, '.', q_strand, '.', utr_attr
					print(*line, sep='\t')
			else:
				print(feature)
		if feature[2]!='CDS':
			print(feature)

# def output_features_minus(q_transcript):#, f):
# 	print(q_transcript)#, file = f)
# 	for feature in qry_db.children(q_transcript, featuretype=['five_prime_UTR','exon','CDS','three_prime_UTR'], order_by='start'): 
# 		print(feature)#, file = f)

def separate_dict(coordinate, frame, seq):
	if frame == 'plus0':
		plus0_dict[coordinate]=seq
	elif frame == 'plus1':
		plus1_dict[coordinate]=seq
	elif frame == 'plus2':
		plus2_dict[coordinate]=seq
	elif frame == 'minus0':
		minus0_dict[coordinate]=seq
	elif frame == 'minus1':
		minus1_dict[coordinate]=seq
	elif frame == 'minus2':
		minus2_dict[coordinate]=seq

plus0_dict={}
plus1_dict={}
plus2_dict={}
minus0_dict={}
minus1_dict={}
minus2_dict={}

#allframe_dict={}

with open(translation_block_gff3, 'r') as in_block:
	for line in in_block:
		if not line.startswith('#'):
			chrom, source, frame, start, end, score, strand, phase, info = line.rstrip().split('\t')[:]
			coordinates_join =chrom, ':', str(start),'-', str(end), strand
			coordinates = ''.join(coordinates_join)
			seq = ''
			if strand == '+':
				seq = fasta.get_seq(chrom, int(start), int(end))
				separate_dict(coordinates, frame, seq)
			else:
				seq = fasta.get_seq(chrom, int(start), int(end), rc=True)
				separate_dict(coordinates, frame, seq)

#for k,v in plus0_dict.items():
#	print (k,v) #V:2092348-2092362+ ACTTCCAAAAGCTAG

print('Parsing..', file = sys.stderr)
plus_ATG = 0
plus_corrected = 0
minus_ATG = 0
minus_corrected = 0

print('##gff-version 3')
for qt_c, q_transcript in enumerate(qry_db.features_of_type('mRNA'),1):
	q_chrom, q_start, q_end, q_strand, q_attributes = q_transcript.seqid, q_transcript.start, q_transcript.end, q_transcript.strand, q_transcript.attributes
	q_transcript_id = q_attributes['Parent'][0]
	genes = qry_db.parents(q_transcript, featuretype='gene', order_by='start')
	coordinates_join =q_chrom, ':', str(q_start),'-', str(q_end)
	coordinates = ''.join(coordinates_join)
	if 'WBGene' in q_transcript_id:
		tracking_q_transcript_id = q_attributes['ID'][0].split(':')[-1] 				#CBG28352.2
		q_transcript_id = tracking_q_transcript_id.split('.')[0]						#CBG28352
		q_gene_id = q_attributes['Parent'][0].split(':')[-1]							#WBGene00089528
		first_cds_seq = ''
		exons = [(CDS.start, CDS.end) for CDS in qry_db.children(q_transcript, featuretype='CDS', order_by='start')]
		if q_transcript.strand != '-':
			(first_cds_s, first_cds_e) = exons[0]
			first_cds_seq=str(fasta.get_seq(q_chrom, first_cds_s, first_cds_e))
			if first_cds_seq.startswith('ATG'):
				plus_ATG+=1
				#print('starts with ATG', tracking_q_transcript_id, coordinates, q_strand, sep = '\t')
				output_features(q_transcript)
				#print(q_transcript, first_cds_seq[:7], 'starts with ATG', tracking_q_transcript_id, coordinates, q_strand, sep = '\t')
			else:
				plus_corrected+=1
				#if tracking_q_transcript_id == 'CBG08558.3':
				corrected_first_cds_seq, corrected_first_seq_s, atg_location = correct_startcodon_plus(q_transcript, first_cds_seq, first_cds_s, first_cds_e)
				print(q_transcript, ' ', first_cds_seq, ' ', corrected_first_cds_seq, '; old first cds coordinate ', q_chrom, ':',first_cds_s, '-', first_cds_e, q_strand, '; new first cds coordinate ', q_chrom, ':',corrected_first_seq_s, '-', first_cds_e, q_strand, sep = '')
				if atg_location == 'upstream':
					output_features_plus_upstream(q_transcript, corrected_first_seq_s)
				if atg_location == 'downstream':
					output_features_plus_downstream(q_transcript, corrected_first_seq_s)
				# print(tracking_q_transcript_id, coordinates, q_strand, first_cds_seq, corrected_first_cds_seq, sep = '\t')
				# print('old first cds coordinate ', q_chrom, ':',first_cds_s, '-', first_cds_e, q_strand, sep = '')
				# print('new first cds coordinate ', q_chrom, ':',corrected_first_seq_s, '-', first_cds_e, q_strand, sep = '')
		else:
			(first_cds_s, first_cds_e) = exons[-1]
			first_cds_seq=str(fasta.get_seq(q_chrom, first_cds_s, first_cds_e, rc=True))
			if first_cds_seq.startswith('ATG'):
				minus_ATG+=1
				#print('starts with ATG', tracking_q_transcript_id, coordinates, q_strand, sep = '\t')
				output_features(q_transcript)
				#print(q_transcript, first_cds_seq[:7], sep = '\t')
				#print(q_transcript, first_cds_seq[:7], 'starts with ATG', tracking_q_transcript_id, coordinates, q_strand, sep = '\t')
			else:
				#corrected_first_cds_seq, corrected_first_seq_e = correct_startcodon_minus(first_cds_seq, first_cds_s, first_cds_e)
				minus_corrected+=1
				corrected_first_cds_seq, corrected_first_seq_e, atg_location =correct_startcodon_minus(q_transcript, first_cds_seq, first_cds_s, first_cds_e)
				print(q_transcript, ' ', first_cds_seq, ' ', corrected_first_cds_seq, '; old first cds coordinate ', q_chrom, ':',first_cds_s, '-', first_cds_e, q_strand, '; new first cds coordinate ', q_chrom, ':',first_cds_s, '-', corrected_first_seq_e, q_strand, sep = '')
				if atg_location == 'upstream':
					output_features_minus_upstream(q_transcript, corrected_first_seq_e)
				if atg_location == 'downstream':
					output_features_minus_downstream(q_transcript, corrected_first_seq_e)
				# print(tracking_q_transcript_id, coordinates, q_strand, first_cds_seq, corrected_first_cds_seq, sep = '\t')
				# print('old first cds coordinate ', q_chrom, ':',first_cds_s, '-', first_cds_e, q_strand, sep = '')
				# print('new first cds coordinate ', q_chrom, ':',first_cds_s, '-', corrected_first_seq_e, q_strand, sep = '')
				# print('\n')
	elif 'XLOC' in q_transcript_id:
		transcript_id = q_attributes['ID'][0]
		exons = [(CDS.start, CDS.end) for CDS in qry_db.children(q_transcript, featuretype='CDS', order_by='start')]
		if q_transcript.strand != '-':
			(first_cds_s, first_cds_e) = exons[0]
			first_cds_seq=str(fasta.get_seq(q_chrom, first_cds_s, first_cds_e))
			if first_cds_seq.startswith('ATG'):
				plus_ATG+=1
				output_features(q_transcript)
			else:
				plus_corrected+=1
				corrected_first_cds_seq, corrected_first_seq_s, atg_location = correct_startcodon_plus(q_transcript, first_cds_seq, first_cds_s, first_cds_e)
				print(q_transcript, ' ', first_cds_seq, ' ', corrected_first_cds_seq, '; old first cds coordinate ', q_chrom, ':',first_cds_s, '-', first_cds_e, q_strand, '; new first cds coordinate ', q_chrom, ':',corrected_first_seq_s, '-', first_cds_e, q_strand, sep = '')
				if atg_location == 'upstream':
					output_features_plus_upstream(q_transcript, corrected_first_seq_s)
				if atg_location == 'downstream':
					output_features_plus_downstream(q_transcript, corrected_first_seq_s)
		else:
			(first_cds_s, first_cds_e) = exons[-1]
			first_cds_seq=str(fasta.get_seq(q_chrom, first_cds_s, first_cds_e, rc=True))
			if first_cds_seq.startswith('ATG'):
				minus_ATG+=1
				output_features(q_transcript)
			else:
				minus_corrected+=1
				corrected_first_cds_seq, corrected_first_seq_e, atg_location =correct_startcodon_minus(q_transcript, first_cds_seq, first_cds_s, first_cds_e)
				print(q_transcript, ' ', first_cds_seq, ' ', corrected_first_cds_seq, '; old first cds coordinate ', q_chrom, ':',first_cds_s, '-', first_cds_e, q_strand, '; new first cds coordinate ', q_chrom, ':',first_cds_s, '-', corrected_first_seq_e, q_strand, sep = '')
				if atg_location == 'upstream':
					output_features_minus_upstream(q_transcript, corrected_first_seq_e)
				if atg_location == 'downstream':
					output_features_minus_downstream(q_transcript, corrected_first_seq_e)

print('plus starts with ATG:', plus_ATG, file = sys.stderr)
print('plus corrected:', plus_corrected, file = sys.stderr)
print('total plus parsed:', plus_ATG+plus_corrected, file = sys.stderr)
print('minus starts with ATG:', minus_ATG, file = sys.stderr)
print('minus corrected:', minus_corrected, file = sys.stderr)
print('total minus parsed:', minus_ATG+minus_corrected, file = sys.stderr)
