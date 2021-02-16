#!/home3/sthio/anaconda/envs/py356/bin/python3.5
#Author: Shinta Thio
#Last updated: Dec 11,2019
# Purpose: 
# To categorize transcripts into 3 general categories and 13 specific categories based on their intron chain.
# 3 general categories: 
# 1) match
# 2) overlap
# 3) novel
# The first two categories are for assembled transcripts that have their intron chains completely match or partially match (i.e., overlapping) those of existing WormBase transcripts. The last category is for assembled transcripts that are not present (i.e., not overlapping) WormBase transcripts.
# 13 specific categories:
# 1) Complete Match (WB confirmed)
# 2) 3' extension
# 3) 5' extension
# 4) 5'&3' extension
# 5) internal within exon
# 6) internal within intron
# 7) alternative donor
# 8) alternative acceptor
# 9) alternative donor&alternative acceptor
# 10) merged
# 11) complex changes
# 12) other:single-exon, partial
# 13) novel

import os, sys, re, pyfaidx
import gffutils as gff
from gffutils import Feature
from collections import defaultdict
from itertools import chain
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.pairwise2 import format_alignment

ref		= sys.argv[1]
qry		= sys.argv[2]
fasta = pyfaidx.Fasta('/home3/sthio/PROJECTS/C_briggsae101/gene_model/c_briggsae.PRJNA10731.WS254.genomic.fa')
match   = sys.argv[3]
r_ext   = sys.argv[4]
l_ext   = sys.argv[5]
lr_ext  = sys.argv[6]
iwe 	= sys.argv[7]
iwi 	= sys.argv[8]
ad  	= sys.argv[9]
aa  	= sys.argv[10]
aaad  	= sys.argv[11]
merged  = sys.argv[12]
complexchanges = sys.argv[13]
novel 	= sys.argv[14]
singleintron = sys.argv[15]
shorter = sys.argv[16]

ref_db_name= ref + '.db'
qry_db_name= qry + '.db'

if os.path.isfile(ref_db_name):
	print('Reference database found. Importing...', file=sys.stderr)
	ref_db = gff.FeatureDB(ref_db_name)
else:
	print('Building reference database...', file=sys.stderr)
	ref_db = gff.create_db(ref, dbfn=ref_db_name,
		merge_strategy='create_unique',
		disable_infer_genes=True,
		disable_infer_transcripts=True,
		verbose=True,
		keep_order=True)

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

def output_features(q_transcript, f):
	parent = qry_db.parents(q_transcript, featuretype='gene')
	for p in parent: 
		print (p, file = f)
	print(q_transcript, file = f)
	for feature in qry_db.children(q_transcript, featuretype=['five_prime_UTR','exon','CDS','three_prime_UTR'], order_by='start'): 
		print(feature, file = f)

def parent_novel(q_transcript, f):
	parent = qry_db.parents(q_transcript, featuretype='gene')
	for p in parent: 
		print (p, file = f)
	output_features(q_transcript, f)

def build_linked_noisoform(tracking_q_transcript_id, temp5):
	for trscrpt, (r_intron_chain, ref_gene_id, ref_coordinates) in ref_transcript_intron_dict.items():
		for i in temp5:
			if i in r_intron_chain: #TO-DO: STRANDEDNESS!!
				if tracking_q_transcript_id not in linked_dict.keys():
					linked_dict[tracking_q_transcript_id]=[coordinates, q_strand, q_transcript_id,  ref_coordinates, trscrpt]
				else:
					v = linked_dict.get(tracking_q_transcript_id)
					if v[-1]!= trscrpt:
						linked_dict[tracking_q_transcript_id].extend((ref_coordinates,trscrpt))
					else:
						pass

def build_linked_isoform(tracking_q_transcript_id, temp5):
	for trscrpt, (r_intron_chain, ref_gene_id, ref_coordinates) in ref_transcript_intron_isoform_dict.items():
		for i in temp5:
			if i in r_intron_chain:
				if tracking_q_transcript_id not in linked_dict_isoform.keys():
					linked_dict_isoform[tracking_q_transcript_id]=[coordinates, q_strand, q_transcript_id, ref_coordinates, trscrpt]
				else:
					v = linked_dict_isoform.get(tracking_q_transcript_id)
					if v[-1]!= trscrpt:
						linked_dict_isoform[tracking_q_transcript_id].extend((ref_coordinates,trscrpt))
					else:
						pass

def gffheader(f):
	header = '##gff-version 3'
	print(header, file=f)

def parse_CDS(q_transcript, id, coordinate, strand, typee):
	cds_chain = ''
	cds_list=[]
	if q_transcript.strand != '-':
		for cds in qry_db.children(q_transcript, featuretype='CDS', order_by='start'):
			cds_seq=Seq(cds.sequence(fasta), generic_dna)
			cds_chain+=cds_seq
	else:
		for cds in qry_db.children(q_transcript, featuretype='CDS', order_by='start'):
			cds_seq=Seq(cds.sequence(fasta), generic_dna)
			cds_list.append(cds_seq)
		for r_cds in reversed(cds_list):
			cds_chain+=r_cds
	cds_transl = cds_chain.translate()
	#print(tracking_q_transcript_id)
	if typee == 'nt': 
		return cds_chain
	elif typee == 'pep':
		return cds_transl

def returnNotMatches(a, b):
	return [[x for x in a if x not in b], [x for x in b if x not in a]]

for g_c, gene in enumerate(ref_db.features_of_type('gene'),1):
	chrom, strand = gene.seqid, gene.strand
	left, right = int(gene.start), int(gene.end)
	#results[(chrom, strand)].append((left, right))

# print(results)
# gene_dict = parse_genes(ref)
# print(gene_dict)
#print(g_c) #21814
#print(len(ref_gene_dict)) #21814

ref_transcript_dict = {}
ref_transcript_intron_dict = {}
ref_transcript_intron_isoform_dict = {}
# ref_transcript_exon_dict = {}
# ref_transcript_exon_isoform_dict = {}
qry_transcript_intron_dict = {}

isoform_c    = 0
no_isoform_c = 0

ref_cbr_pep_dict={}
ref_cel_pep_dict={}
ref_cbr_pep_c = 0  
ref_cel_pep_c = 0

cbr_cel_ortho_dict={}
cel_cbr_ortho_dict={}

ref_cel_pep_34_dict = {}

# with open(cbr_cel_ortholog_tsv_in, 'r') as ortho_info:
# 	for line in ortho_info:
# 		WBGene, CBG, cbr_alias, biotype, cel_ortholog = line.rstrip().split('\t')[:]
# 		cbr_cel_ortho_dict[CBG] = [WBGene, CBG, cbr_alias, biotype, cel_ortholog]
# 		if cbr_alias.startswith('Cbr'):
# 			cel_alias=cbr_alias[4:]
# 			cel_cbr_ortho_dict[cel_alias]=cbr_alias
# 		else:
# 			cel_alias=cel_ortholog
# 			cel_cbr_ortho_dict[cel_alias]=cbr_alias

# count = 0
# count_none = 0

# PID_dict = {}
# with open(cbr_ref_pep_in, 'r') as in_cbr_ref, open(cel_ref_pep_in, 'r') as in_cel_ref:
# 	#parse_peptide(in_cbr_ref, ref_cbr_pep_dict)
# 	#parse_peptide(in_cel_ref, ref_cel_pep_dict)
# 	for record in SeqIO.parse(in_cbr_ref, 'fasta'):
# 		seqid, pep = record.id, record.seq
# 		cel_ortholog = ''
# 		if 'Cbr' in record.description: #6,635
# 			locus=record.description.split('\t')[3]
# 			cel_ortholog=locus.split('locus:Cbr-')[1] #print(locus, cel_ortholog)
# 		else:
# 			cel_ortholog='other' #None or hypothetical
# 		ref_cbr_pep_c+=1
# 		ref_cbr_pep_dict[seqid]=pep, cel_ortholog if seqid not in ref_cbr_pep_dict else print ("Warning: redundant transcript found.")
# 	for record in SeqIO.parse(in_cel_ref, 'fasta'):
# 		seqid, pep = record.id, record.seq
# 		cel_ortholog = ''
# 		if 'locus' in record.description: #14,836
# 			locus=record.description.split('\t')[3]
# 			cel_ortholog=locus.split('locus:')[1] #print(locus, cel_ortholog)
# 		else:
# 			cel_ortholog='other'
# 		ref_cel_pep_c+=1
# 		ref_cel_pep_dict[seqid]=pep, cel_ortholog if seqid not in ref_cel_pep_dict else print ("Warning: redundant transcript found.")
# 		ref_cel_pep_34_dict[cel_ortholog]=pep if cel_ortholog != 'other' else 0
# 	for cbr_transcript, (cbr_pep, cel_ortho1) in ref_cbr_pep_dict.items():
# 		if cel_ortho1 != 'other':
# 			cel_pep = ref_cel_pep_34_dict.get(cel_ortho1)
# 			if cel_pep != None:
# 				count+=1
# 				# l_alignments = pairwise2.align.localxx(cbr_pep, cel_pep)
# 				# len_ref_match = 0
# 				# for a in l_alignments:
# 				# 	len_ref_match = int(a[2])
# 				# ori_PID = len_ref_match/len(cel_pep)*100
# 				# PID_dict[cbr_transcript]=cel_ortho1, ori_PID
# 			else:
# 				count_none+=1

# for cbr_CBG, (cel_ortho1, oriPID) in PID_dict.items():
# 	print(cbr_CBG, cel_ortho1, oriPID,sep='\t') #done - file: cbr_cel_original_PID.tsv
# ori_PID_dict={}

# with open(original_PID) as ori:
# 	for line in ori:
# 		cbr_transcript, cel_ortholog, ori_PID = line.rstrip().split('\t')[:]
# 		ori_PID_dict[cbr_transcript]=cel_ortholog, ori_PID

def calc_rev_PID(q_transcript_id, coordinates, strand):
	if q_transcript_id in ori_PID_dict.keys():
		cel_ortholog = ref_cbr_pep_dict.get(q_transcript_id)[1]
		cel_pep = ref_cel_pep_34_dict.get(cel_ortholog)
		rev_cbr_pep = parse_CDS(q_transcript, tracking_q_transcript_id, coordinates, q_strand, 'pep')
		l_alignments = pairwise2.align.localxx(rev_cbr_pep, cel_pep)
		len_ref_match = 0
		for a in l_alignments:
			#print(format_alignment(*a))
			len_ref_match = int(a[2])
		#print('cbr', rev_cbr_pep, len(rev_cbr_pep))
		#print('cel', cel_pep, len(cel_pep))
		#print('-------------------------', len_ref_match)
		ori_PID = float(ori_PID_dict.get(q_transcript_id)[1])
		rev_PID = len_ref_match/len(cel_pep)*100
		delta_PID = rev_PID - ori_PID
		return ori_PID, rev_PID, delta_PID
		#print(category, coordinates, strand, 'ori PID', "%.2f" %ori_PID, 'rev PID', "%.2f" %rev_PID, 'delta_PID', "%.2f" %delta_PID)
	else:
		return '-','-','-'

#print(count, count_none)
#print(ref_cbr_pep_c, ref_cel_pep_c) #21863 28071

''' parsing the c. briggsae protein coding transcripts '''
for t_c, transcript in enumerate(ref_db.features_of_type('mRNA'),1):
	m_chrom, m_start, m_end, m_strand, m_attributes = transcript.seqid, transcript.start, transcript.end, transcript.strand, transcript.attributes
	transcript_id = transcript.attributes['ID'][0].split(':')[-1]
	ref_transcript_dict[transcript_id]=(m_start,m_end,m_strand)
	gene_id = transcript.attributes['Parent'][0].split(':')[-1]
	coordinates_join =m_chrom, ':', str(m_start),'-', str(m_end), m_strand
	coordinates = ''.join(coordinates_join)
	#exon_chain=[(exon.start, exon.end) for exon in ref_db.children(transcript, featuretype='CDS', order_by='start')]
	#intron_chain = [(intron.start, intron.end) for intron in ref_db.children(transcript, featuretype='intron', order_by='start')]
	exons = [(CDS.start, CDS.end) for CDS in ref_db.children(transcript, featuretype='CDS', order_by='start')]
	intron_chain = [(exons[n][1]+1, exons[n+1][0]-1) for n in range(len(exons)-1)]
	''' separating if those protein-coding transcripts have single or multiple isoforms per gene''' 
	''' with multiple isoforms'''
	if '.' in transcript_id:
		isoform_exons = [(CDS.start, CDS.end) for CDS in ref_db.children(transcript, featuretype='CDS', order_by='start')]
		set_isoform_exons = set(isoform_exons)
		list_set_isoform_exons = list(set_isoform_exons)
		list_set_isoform_exons.sort()
		isoform_intron_chain = [(list_set_isoform_exons[n][1]+1, list_set_isoform_exons[n+1][0]-1) for n in range(len(list_set_isoform_exons)-1)]
		ref_transcript_intron_isoform_dict[transcript_id]=isoform_intron_chain, gene_id, coordinates      ##1035
		isoform_c+=1	#isoforms: 1035
	''' with single isoform per gene '''
	else:
		ref_transcript_intron_dict[transcript_id]=intron_chain, gene_id, coordinates				
		no_isoform_c+=1#no isoforms: 21,440

#print(ref_transcript_dict)
#print(len(ref_transcript_dict))
#print(len(ref_transcript_intron_isoform_dict)) 
#total protein-coding genes = 21,814; total mRNAs = 22,745
#print(isoform_c, no_isoform_c) #1035 21440 - total all=22,745
# print(t_c)										#22,475
# print(len(ref_transcript_dict))					#22,475
# print(len(ref_transcript_intron_dict))			#21,440
# print(len(ref_transcript_intron_isoform_dict))	 #1,035

q_notxloc_c= 0
q_xloc_c	= 0
q_isoform_c= 0
q_noisoform_c= 0

ic_same_c	= 0
ic_longer_c= 0
ic_shorter_c= 0

ic_same_completelymatch_c = 0
ic_same_noCDStranscript_c = 0
ic_same_differentcontent_c = 0
ic_same_diff_singleintron_c = 0
ic_same_diff_singleintron_AA_c = 0
ic_same_diff_singleintron_AD_c = 0
ic_same_diff_singleintron_both_AA_AD_c = 0
ic_same_diff_notsingleintron_c = 0
ic_same_diff_notsingleintron_1diff_c = 0
ic_same_diff_notsingleintron_1diff_AA_c = 0
ic_same_diff_notsingleintron_1diff_AD_c = 0
ic_same_diff_notsingleintron_1diff_both_AA_AD_c = 0
ic_same_diff_notsingleintron_morediff_c = 0

ic_longer_ref_subset_qry_c = 0
ic_longer_ref_subset_qry_refiszero_c = 0
ic_longer_ref_subset_qry_refnotzero_c = 0
ic_longer_ref_subset_qry_refnotzero_left_extension_c = 0
ic_longer_ref_subset_qry_refnotzero_right_extension_c = 0
ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c = 0
ic_longer_ref_subset_qry_refnotzero_complexchanges_ext_iwe_c =0
ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c = 0
ic_longer_ref_subset_qry_refnotzero_other_c = 0
ic_longer_ref_subset_qry_other_c = 0
ic_longer_ref_notsubset_qry_c = 0
ic_longer_ref_notsubset_qry_internalwithinintron_c = 0
ic_longer_ref_notsubset_qry_complexchanges_c = 0
ic_longer_ref_notsubset_qry_other_c = 0

ic_longer_ref_subset_qry_refnotzero_right_linked_c = 0
ic_longer_ref_subset_qry_refnotzero_left_linked_c = 0
ic_longer_ref_subset_qry_refnotzero_left_right_linked_c = 0
ic_longer_ref_subset_qry_other_linked_c = 0

i_ic_same_c= 0
i_ic_longer_c= 0
i_ic_shorter_c= 0

i_ic_same_noCDStranscript_c = 0
i_ic_same_completelymatch_c = 0
i_ic_same_differentcontent_c = 0
i_ic_same_diff_singleintron_c = 0
i_ic_same_diff_singleintron_AD_c = 0
i_ic_same_diff_singleintron_AA_c = 0 
i_ic_same_diff_singleintron_both_AA_AD_c = 0
i_ic_same_diff_notsingleintron_c = 0
i_ic_same_diff_notsingleintron_1diff_c = 0
i_ic_same_diff_notsingleintron_1diff_AA_c = 0
i_ic_same_diff_notsingleintron_1diff_AD_c = 0
i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c = 0
i_ic_same_diff_notsingleintron_morediff_c = 0

i_ic_longer_ref_subset_qry_c = 0
i_ic_longer_ref_subset_qry_refiszero_c = 0
i_ic_longer_ref_subset_qry_refnotzero_c = 0
i_ic_longer_ref_subset_qry_refnotzero_left_extension_c = 0
i_ic_longer_ref_subset_qry_refnotzero_right_extension_c = 0
i_ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c = 0
i_ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c = 0
i_ic_longer_ref_subset_qry_refnotzero_other_c = 0
i_ic_longer_ref_subset_qry_other_c = 0
i_ic_longer_ref_notsubset_qry_c = 0
i_ic_longer_ref_notsubset_qry_internalwithinintron_c = 0
i_ic_longer_ref_notsubset_qry_complexchanges_c = 0
i_ic_longer_ref_notsubset_qry_other_c = 0

i_ic_longer_ref_subset_qry_refnotzero_left_linked_c = 0
i_ic_longer_ref_subset_qry_refnotzero_right_linked_c = 0
#i_ic_longer_ref_subset_qry_refnotzero_left_right_linked_c = 0
i_ic_longer_ref_subset_qry_other_linked_c = 0

linked_dict = {}
linked_dict_isoform = {}
isoform_dict = {}
isoform_zoomout_dict = {}

c_conserved = 0
with open(match, 'w') as o_match, open(r_ext, 'w') as o_r_ext, open(l_ext, 'w') as o_l_ext, open(lr_ext, 'w') as o_lr_ext, open(iwe, 'w') as o_iwe, open(iwi, 'w') as o_iwi, open(ad, 'w') as o_ad, open(aa, 'w') as o_aa, open(aaad, 'w') as o_aaad, open(merged, 'w') as o_merged, open(complexchanges, 'w') as o_cc, open(novel, 'w') as o_novel, open(singleintron, 'w') as o_si, open(shorter, 'w') as o_short:
	o_files = [o_match, o_r_ext, o_l_ext, o_lr_ext, o_iwe, o_iwi, o_ad, o_aa, o_aaad, o_merged, o_cc, o_novel, o_si, o_short]
	for o_file in o_files:
		gffheader(o_file)
#with open(overlap_and_novel, 'w') as o_overlap_novel, open(overlap_only, 'w') as o_overlap_only, open(overlap_nomerged, 'w') as o_overlap_nomerged, open(merged_only, 'w') as o_merged:
#with open(novel_only,'w') as o_novel_only:
#with open(r_ext,'w') as r_ext_pid, open(l_ext, 'w') as l_ext_pid, open(lr_ext, 'w') as lr_ext_pid, open(iwe, 'w') as iwe_pid, open(iwi, 'w') as iwi_pid, open(ad, 'w') as ad_pid, open(aa, 'w') as aa_pid, open(aaad, 'w') as aaad_pid, open(complexchanges, 'w') as complexchanges_pid:
	''' parsing the RNA-seq protein coding transcripts '''
	for qt_c, q_transcript in enumerate(qry_db.features_of_type('mRNA'),1):
		q_chrom, q_start, q_end, q_strand, q_attributes = q_transcript.seqid, q_transcript.start, q_transcript.end, q_transcript.strand, q_transcript.attributes
		''' transdecoder trasncript format: 
		III     transdecoder    mRNA    5929166 5935467 .       -       .       Parent=Gene:WBGene00033997;Name=ORF type:complete len:780 (+)%2Cscore%3D167.38;ID=Transcript:CBG13202.3.p1 '''
		q_transcript_id = q_attributes['Parent'][0]
		genes = qry_db.parents(q_transcript, featuretype='gene', order_by='start')
		coordinates_join =q_chrom, ':', str(q_start),'-', str(q_end)
		coordinates = ''.join(coordinates_join)
		if 'WBGene' in q_transcript_id:
			q_notxloc_c +=1
			tracking_q_transcript_id = q_attributes['ID'][0].split(':')[-1] 				#CBG28352.2
			q_transcript_id = tracking_q_transcript_id.split('.')[0]						#CBG28352
			q_gene_id = q_attributes['Parent'][0].split(':')[-1]							#WBGene00089528
			exons = [(CDS.start, CDS.end) for CDS in qry_db.children(q_transcript, featuretype='CDS', order_by='start')]
			intron_chain = [(exons[n][1]+1, exons[n+1][0]-1) for n in range(len(exons)-1)]
			qry_transcript_intron_dict[tracking_q_transcript_id]=intron_chain
			if q_transcript_id in ref_transcript_intron_dict.keys(): 						##23975
				q_noisoform_c+=1
				ref_transcript_id = q_transcript_id
				ref_intron_chain = ref_transcript_intron_dict.get(q_transcript_id)[0]
				set_intron_chain = set(intron_chain)
				set_ref_intron_chain = set(ref_intron_chain)
				if len(intron_chain) == len(ref_intron_chain):
					ic_same_c+=1															#10,312 --- 10,665
					if len(intron_chain)==0 and len(ref_intron_chain)==0:
						ic_same_noCDStranscript_c+=1										#    19 ---     21
						category = "Ref CDS intron absent"
						output_features(q_transcript, o_si)
						#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category, sep= '\t', file = o_si)
					elif intron_chain == ref_intron_chain:
						ic_same_completelymatch_c+=1										# 7,523 ---  7,825 #COMPLETELY MATCH 
						category = "Complete match"
						#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
						# #print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t')
						#if delta_PID == 0.0:
							#c_conserved+=1 #2709
						# if 'CBG06882' in q_transcript_id:
						# 	parse_CDS(q_transcript, tracking_q_transcript_id, coordinates, q_strand, 'nt')
						# 	parse_CDS(q_transcript, tracking_q_transcript_id, coordinates, q_strand, 'pep')
						output_features(q_transcript, o_match)
					else:
						ic_same_differentcontent_c+=1											# 2,770 ---  2,819 
						if len(intron_chain) == 1:	
							ic_same_diff_singleintron_c+=1									#   89 (AA:30, AD:32, both:27) --- 93 (AA:30, AD:34, both:29)
							if intron_chain[0][0] == ref_intron_chain[0][0]:
								if q_strand == '-':
									category = "Alt. donor"
									ic_same_diff_singleintron_AD_c+=1
									#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
									#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
									output_features(q_transcript, o_ad)
								else:
									category = "Alt. acceptor"
									ic_same_diff_singleintron_AA_c+=1
									#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
									#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
									output_features(q_transcript, o_aa)
							elif intron_chain[0][-1] == ref_intron_chain[0][-1]:
								if q_strand == '-':
									category = "Alt. acceptor"
									ic_same_diff_singleintron_AA_c+=1
									#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
									#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
									output_features(q_transcript, o_aa)
								else:
									category = "Alt. donor"
									ic_same_diff_singleintron_AD_c+=1
									#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
									#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
									output_features(q_transcript, o_ad)
							else:
								category = "Alt. acceptor and donor"
								ic_same_diff_singleintron_both_AA_AD_c +=1
								#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
								#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aaad_pid)
								output_features(q_transcript, o_aaad)
						else:
							ic_same_diff_notsingleintron_c+=1								# 2,681 ---  2,726
							temp3 = sorted([x for x in set_intron_chain if x not in set_ref_intron_chain])
							temp4 = sorted([y for y in set_ref_intron_chain if y not in set_intron_chain])
							if len(temp3) ==1:
								ic_same_diff_notsingleintron_1diff_c+=1					#1,845 (AA:816, AD:728, both:301)
								if temp3[0][0]==temp4[0][0] and temp3[0][-1]!=temp4[0][-1]:
									if q_strand == '-':
										category = "Alt. donor"
										ic_same_diff_notsingleintron_1diff_AD_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
										output_features(q_transcript, o_ad)
									else:
										category = "Alt. acceptor"
										ic_same_diff_notsingleintron_1diff_AA_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
										output_features(q_transcript, o_aa)
								elif temp3[0][0]!=temp4[0][0] and temp3[0][-1]==temp4[0][-1]:
									if q_strand == '-':
										category = "Alt. acceptor"
										ic_same_diff_notsingleintron_1diff_AA_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
										output_features(q_transcript, o_aa)
									else:
										category = "Alt. donor"
										ic_same_diff_notsingleintron_1diff_AD_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
										output_features(q_transcript, o_ad)
								else:
									category = "Alt. acceptor and donor"
									ic_same_diff_notsingleintron_1diff_both_AA_AD_c+=1											 # 1,845 --- 1,885
									#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
									#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aaad_pid)
									output_features(q_transcript, o_aaad)
							else:														
								category = "Complex changes, multiple AS events" 
								#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_overlap_nomerged)
								#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category, sep= '\t', file = tbd_out) #fileout: TBD_existing_transcripts.tsv
								ic_same_diff_notsingleintron_morediff_c+=1
													 #   836 --- 841
								#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
								#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = complexchanges_pid)
								output_features(q_transcript, o_cc)
				elif len(intron_chain)>len(ref_intron_chain):								 # 2,908 --- 3,215
					ic_longer_c+=1
					temp5 = sorted([x for x in set_intron_chain if x not in set_ref_intron_chain])
					temp6 = sorted([y for y in set_ref_intron_chain if y not in set_intron_chain])
					if set_ref_intron_chain.issubset(set_intron_chain):					 # 1,225 --- 1,482
						ic_longer_ref_subset_qry_c+=1			
						if len(ref_intron_chain)==0:
							category = 'Ref CDS intron absent, qry present'
							ic_longer_ref_subset_qry_refiszero_c+=1								#     74 ---    88
							output_features(q_transcript, o_si)
							#output_features(q_transcript)
							#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category, sep= '\t', file = o_si)
						elif len(set_ref_intron_chain) > 0:
							ic_longer_ref_subset_qry_refnotzero_c+=1						#  1,151 (total ext: 800; 5' ext:458, 3' ext:310, both ext:32) --- 1,394 (total ext: 1,039; 5' ext: 688, 3' ext: 313, bth ext: 38)
							start_pos_intron_chain = intron_chain[0]
							end_pos_intron_chain = intron_chain[-1] 
							start_pos_ref_intron_chain = ref_intron_chain[0]
							end_pos_ref_intron_chain= ref_intron_chain[-1] 
							if start_pos_intron_chain == start_pos_ref_intron_chain and end_pos_intron_chain != end_pos_ref_intron_chain:
								if q_strand == '-':
									build_linked_noisoform(tracking_q_transcript_id, temp5)
									if tracking_q_transcript_id in linked_dict.keys():
										category = "Linked, 5'"
										ic_longer_ref_subset_qry_refnotzero_left_linked_c+=1
										k, v = tracking_q_transcript_id, linked_dict.get(tracking_q_transcript_id)
										output_features(q_transcript, o_merged)
										#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
										# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
										#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
									else:
										category = "5' extension"
										ic_longer_ref_subset_qry_refnotzero_left_extension_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = l_ext_pid)
										output_features(q_transcript, o_l_ext)
								else:
									build_linked_noisoform(tracking_q_transcript_id, temp5)
									if tracking_q_transcript_id in linked_dict.keys():
										category = "Linked, 3'"
										ic_longer_ref_subset_qry_refnotzero_right_linked_c+=1
										output_features(q_transcript, o_merged)
										k, v = tracking_q_transcript_id, linked_dict.get(tracking_q_transcript_id)
										#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_merged)
										#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
										# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
										#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
									else:
										category = "3' extension"
										ic_longer_ref_subset_qry_refnotzero_right_extension_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = r_ext_pid)
										output_features(q_transcript, o_r_ext)
							elif end_pos_intron_chain == end_pos_ref_intron_chain and start_pos_intron_chain != start_pos_ref_intron_chain:
								if q_strand == '-':
									build_linked_noisoform(tracking_q_transcript_id, temp5)
									if tracking_q_transcript_id in linked_dict.keys():
										category = "Linked, 3'"
										ic_longer_ref_subset_qry_refnotzero_right_linked_c+=1
										output_features(q_transcript, o_merged)
										k, v = tracking_q_transcript_id, linked_dict.get(tracking_q_transcript_id)
										#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_merged)
										#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
										# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
										#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
									else:
										category = "3' extension"
										ic_longer_ref_subset_qry_refnotzero_right_extension_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = r_ext_pid)
										output_features(q_transcript, o_r_ext)
								else:
									build_linked_noisoform(tracking_q_transcript_id, temp5)
									if tracking_q_transcript_id in linked_dict.keys():
										category = "Linked, 5'"
										ic_longer_ref_subset_qry_refnotzero_left_linked_c+=1
										output_features(q_transcript, o_merged)
										k, v = tracking_q_transcript_id, linked_dict.get(tracking_q_transcript_id)
										#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_merged)
										#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
										# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
										#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
									else:
										category = "5' extension"
										ic_longer_ref_subset_qry_refnotzero_left_extension_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = l_ext_pid)
										output_features(q_transcript, o_l_ext)
							elif start_pos_intron_chain != start_pos_ref_intron_chain and end_pos_intron_chain != end_pos_ref_intron_chain:
								build_linked_noisoform(tracking_q_transcript_id, temp5)
								if tracking_q_transcript_id in linked_dict.keys():
									category = "linked, additional intron on 5' & 3'"
									ic_longer_ref_subset_qry_refnotzero_left_right_linked_c+=1
									k, v = tracking_q_transcript_id, linked_dict.get(tracking_q_transcript_id)
									output_features(q_transcript, o_merged)
									#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_merged)
									#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
									# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
									#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
								else:
									if len(intron_chain) == len(ref_intron_chain)+2:
										category = "5' and 3' extension"
										ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = lr_ext_pid)
										output_features(q_transcript, o_lr_ext)
									else:
										category = "complex changes (both extension+iwe)"
										ic_longer_ref_subset_qry_refnotzero_complexchanges_ext_iwe_c+=1
										#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
										#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = complexchanges_pid)
										output_features(q_transcript, o_cc)
							elif start_pos_intron_chain == start_pos_ref_intron_chain and end_pos_intron_chain == end_pos_ref_intron_chain:
								category = "Internal within exon"
								ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c+=1      #  351 ---   355
								#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
								#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = iwe_pid)
								output_features(q_transcript, o_iwe)
								# for i in temp5:
								# 	delta = (i[1]-i[0]+1)
								# 	frame = delta/3
								# 	if delta%3==0:
								# 		continue #print('delta=', delta, 'div', frame)
								# 	else:
								# 		print('delta=', delta, 'div', frame, q_gene_id, tracking_q_transcript_id, coordinates, q_strand)
							else:
								category = "other"
								ic_longer_ref_subset_qry_refnotzero_other_c +=1			 #    0 ---      0
						else:
							ic_longer_ref_subset_qry_other_c+=1							 #    0 ---      0
					else:
						ic_longer_ref_notsubset_qry_c+=1							#1,683
						start_pos_intron_chain = intron_chain[0]
						s_start_pos_intron_chain = intron_chain[0][0]
						end_pos_intron_chain = intron_chain[-1] 
						e_end_pos_intron_chain = intron_chain[-1][-1] 
						start_pos_ref_intron_chain = ref_intron_chain[0]
						s_start_pos_ref_intron_chain = ref_intron_chain[0][0]
						e_end_pos_ref_intron_chain= ref_intron_chain[-1][-1]
						if s_start_pos_intron_chain ==s_start_pos_ref_intron_chain and e_end_pos_intron_chain == e_end_pos_ref_intron_chain: #CHANGES ARE ALWAYS INTERNAL OF LEFTMOST AND RIGHTMOST OF THE TRANSCRIPT!
							if temp5[0][0] == temp6[0][0] and temp5[-1][-1] ==temp6[0][-1]:
								category = "Internal transcript boundaries, within intron" # additional intron(S) within intron=96						
								ic_longer_ref_notsubset_qry_internalwithinintron_c+=1 #NOT LIMITED TO ONE INTRON ONLY
								#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
								#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = iwi_pid)
								output_features(q_transcript, o_iwi)
							else:
								category = "Complex changes, internal transcript boundaries" #internal and/or other modifications
								ic_longer_ref_notsubset_qry_complexchanges_c+=1	# 564
								#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
								#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = complexchanges_pid)
								output_features(q_transcript, o_cc)
						else:
							build_linked_noisoform(tracking_q_transcript_id, temp5)
							if tracking_q_transcript_id in linked_dict.keys():
								category = "linked"
								ic_longer_ref_subset_qry_other_linked_c+=1
								output_features(q_transcript, o_merged)
								k, v = tracking_q_transcript_id, linked_dict.get(tracking_q_transcript_id)
								#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_merged)
								#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
								# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
								#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
							else:
								category = "Complex changes, diff transcript boundaries" 
								ic_longer_ref_notsubset_qry_other_c+=1					#1,027
								#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
								#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = complexchanges_pid)
								output_features(q_transcript, o_cc)
				elif len(intron_chain)< len(ref_intron_chain):
					category = 'Shorter'
					ic_shorter_c+=1											# 10,755
					output_features(q_transcript, o_short)
			else:
				for k,(ref_intron_chain, ref_gene_id, ref_coordinates) in ref_transcript_intron_isoform_dict.items():
					k_modif=k.split('.')[0]
					if q_transcript_id == k_modif:
						q_isoform_c+=1							 ## 1355 (lotsa doubles since each assembled transcript isoforms are being compared to each of the reference transcript isoforms)
						ref_transcript_id = q_transcript_id
						set_intron_chain = set(intron_chain)
						set_ref_intron_chain = set(ref_intron_chain)
						if len(intron_chain) == len(ref_intron_chain):
							i_ic_same_c+=1
							if len(intron_chain)==0 and len(ref_intron_chain)==0:
								if tracking_q_transcript_id not in isoform_dict.keys():
									isoform_dict[tracking_q_transcript_id]=1
									i_ic_same_noCDStranscript_c+=1			#22 with doubles (wd)
									category = "Ref CDS intron absent"
									output_features(q_transcript, o_si)
									#output_features(q_transcript)
									#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category, sep= '\t', file = o_si)
							elif intron_chain == ref_intron_chain:
								#i_ic_same_completelymatch_c+=1 	        #621 (wd)
								category = "Complete match"
								if tracking_q_transcript_id not in isoform_dict.keys():
									isoform_dict[tracking_q_transcript_id]=1
									i_ic_same_completelymatch_c+=1
									##print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t')
									#if delta_PID == 0.0:
										#c_conserved+=1 #2811-2709
									output_features(q_transcript, o_match)
							else:
								i_ic_same_differentcontent_c+=1	 	   
								if len(intron_chain) == 1:
									i_ic_same_diff_singleintron_c+=1    #2 (both AA&AD:2)
									if intron_chain[0][0] == ref_intron_chain[0][0]:
										if q_strand == '-':
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												category = "Alt. donor"
												i_ic_same_diff_singleintron_AD_c+=1
												#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
												#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
												output_features(q_transcript, o_ad)
										else:
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												category = "Alt. acceptor"
												i_ic_same_diff_singleintron_AA_c+=1
												#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
												#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
												output_features(q_transcript, o_aa)
									elif intron_chain[0][-1] == ref_intron_chain[0][-1]:
										if q_strand == '-':
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												category = "Alt. acceptor"
												i_ic_same_diff_singleintron_AA_c+=1
												#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
												#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
												output_features(q_transcript, o_aa)
										else:
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												category = "Alt. donor"
												i_ic_same_diff_singleintron_AD_c+=1
												#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
												#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
												output_features(q_transcript, o_ad)
									else:
										if tracking_q_transcript_id not in isoform_dict.keys():
											isoform_dict[tracking_q_transcript_id]=1
											category = "Alt. acceptor and donor"
											i_ic_same_diff_singleintron_both_AA_AD_c +=1
											#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
											#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aaad_pid)
											output_features(q_transcript, o_aaad)
								else:
									i_ic_same_diff_notsingleintron_c+=1  
									temp3 = [x for x in set_intron_chain if x not in set_ref_intron_chain]
									temp4 = [y for y in set_ref_intron_chain if y not in set_intron_chain]
									if len(temp3) ==1:
										i_ic_same_diff_notsingleintron_1diff_c+=1					#58 (AA:42, AD:12, both AA&AD:4)
										if temp3[0][0]==temp4[0][0] and temp3[0][-1]!=temp4[0][-1]:
											if q_strand == '-':
												if tracking_q_transcript_id not in isoform_dict.keys():
													isoform_dict[tracking_q_transcript_id]=1
													category = "Alt. donor"
													i_ic_same_diff_notsingleintron_1diff_AD_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
													output_features(q_transcript, o_ad)
											else:
												if tracking_q_transcript_id not in isoform_dict.keys():
													isoform_dict[tracking_q_transcript_id]=1
													category = "Alt. acceptor"
													i_ic_same_diff_notsingleintron_1diff_AA_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
													output_features(q_transcript, o_aa)
										elif temp3[0][0]!=temp4[0][0] and temp3[0][-1]==temp4[0][-1]:
											if q_strand == '-':
												if tracking_q_transcript_id not in isoform_dict.keys():
													isoform_dict[tracking_q_transcript_id]=1
													category = "Alt. acceptor"
													i_ic_same_diff_notsingleintron_1diff_AA_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aa_pid)
													output_features(q_transcript, o_aa)
											else:
												if tracking_q_transcript_id not in isoform_dict.keys():
													isoform_dict[tracking_q_transcript_id]=1
													category = "Alt. donor"
													i_ic_same_diff_notsingleintron_1diff_AD_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = ad_pid)
													output_features(q_transcript, o_ad)
										else:
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												category = "Alt. acceptor and donor"
												i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c+=1
												#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
												#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = aaad_pid)
												output_features(q_transcript, o_aaad)
									else:
										if tracking_q_transcript_id not in isoform_dict.keys():
											isoform_dict[tracking_q_transcript_id]=1											
											category = "Complex changes, multiple AS events" #same length, diff content, complex changes' 
											i_ic_same_diff_notsingleintron_morediff_c+=1			#18
											#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
											#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = complexchanges_pid)
											output_features(q_transcript, o_cc)
						elif len(intron_chain)>len(ref_intron_chain):
							i_ic_longer_c+=1
							temp5 = sorted([x for x in set_intron_chain if x not in set_ref_intron_chain])
							temp6 = sorted([y for y in set_ref_intron_chain if y not in set_intron_chain])
							if set_ref_intron_chain.issubset(set_intron_chain):			#1,225
								i_ic_longer_ref_subset_qry_c+=1			
								if len(ref_intron_chain)==0:
									if tracking_q_transcript_id not in isoform_dict.keys():
										isoform_dict[tracking_q_transcript_id]=1
										category = 'Ref CDS intron absent, qry present'
										i_ic_longer_ref_subset_qry_refiszero_c+=1						#   74
										output_features(q_transcript, o_si)
								elif len(set_ref_intron_chain) > 0:
									i_ic_longer_ref_subset_qry_refnotzero_c+=1				#1,151
									start_pos_intron_chain = intron_chain[0]
									end_pos_intron_chain = intron_chain[-1] 
									start_pos_ref_intron_chain = ref_intron_chain[0]
									end_pos_ref_intron_chain= ref_intron_chain[-1] 
									if start_pos_intron_chain == start_pos_ref_intron_chain and end_pos_intron_chain != end_pos_ref_intron_chain:
										if q_strand == '-':
											build_linked_isoform(tracking_q_transcript_id, temp5)
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												if tracking_q_transcript_id in linked_dict_isoform.keys():
													category = "Linked, 5'"
													i_ic_longer_ref_subset_qry_refnotzero_left_linked_c+=1
													output_features(q_transcript, o_merged)
													k, v = tracking_q_transcript_id, linked_dict_isoform.get(tracking_q_transcript_id)
													#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_merged)
													#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
													# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
													#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
												else:
													category = "5' extension"
													i_ic_longer_ref_subset_qry_refnotzero_left_extension_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = l_ext_pid)
													output_features(q_transcript, o_l_ext)
										else:
											build_linked_isoform(tracking_q_transcript_id, temp5)
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												if tracking_q_transcript_id in linked_dict_isoform.keys():
													category = "Linked, 3'"
													i_ic_longer_ref_subset_qry_refnotzero_right_linked_c+=1
													output_features(q_transcript, o_merged)
													k, v = tracking_q_transcript_id, linked_dict_isoform.get(tracking_q_transcript_id)
													#--output_features(q_transcript, o_overlap_novel), output_features(q_transcript, o_overlap_only), output_features(q_transcript, o_merged)
													#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
													# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
													#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
												else:
													category = "3' extension"
													i_ic_longer_ref_subset_qry_refnotzero_right_extension_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = r_ext_pid)
													output_features(q_transcript, o_r_ext)
									elif end_pos_intron_chain == end_pos_ref_intron_chain and start_pos_intron_chain != start_pos_ref_intron_chain:
										if q_strand == '-':
											build_linked_isoform(tracking_q_transcript_id, temp5)
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												if tracking_q_transcript_id in linked_dict_isoform.keys():
													category = "Linked, 3'"
													i_ic_longer_ref_subset_qry_refnotzero_right_linked_c+=1
													output_features(q_transcript, o_merged)
													k, v = tracking_q_transcript_id, linked_dict_isoform.get(tracking_q_transcript_id)
													#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
													# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
													#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
												else:
													category = "3' extension"
													i_ic_longer_ref_subset_qry_refnotzero_right_extension_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = r_ext_pid)
													output_features(q_transcript, o_r_ext)
										else:
											build_linked_isoform(tracking_q_transcript_id, temp5)
											if tracking_q_transcript_id not in isoform_dict.keys():
												isoform_dict[tracking_q_transcript_id]=1
												if tracking_q_transcript_id in linked_dict_isoform.keys():
													category = "Linked, 5'"
													i_ic_longer_ref_subset_qry_refnotzero_left_linked_c+=1
													output_features(q_transcript, o_merged)
													k, v = tracking_q_transcript_id, linked_dict_isoform.get(tracking_q_transcript_id)
													#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
													# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
													#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
												else:
													category = "5' extension"
													i_ic_longer_ref_subset_qry_refnotzero_left_extension_c+=1
													#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
													#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = l_ext_pid)
													output_features(q_transcript, o_l_ext)
									elif start_pos_intron_chain != start_pos_ref_intron_chain and end_pos_intron_chain != end_pos_ref_intron_chain:
										if tracking_q_transcript_id not in isoform_dict.keys():
											isoform_dict[tracking_q_transcript_id]=1
											category = "5' and 3' extension"
											i_ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c+=1
											#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
											#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = lr_ext_pid)
											output_features(q_transcript, o_lr_ext)
									elif start_pos_intron_chain == start_pos_ref_intron_chain and end_pos_intron_chain == end_pos_ref_intron_chain:
										if tracking_q_transcript_id not in isoform_dict.keys():
											isoform_dict[tracking_q_transcript_id]=1
											category = "Internal within exon"
											i_ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c+=1
											#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
											#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = iwe_pid)
											output_features(q_transcript, o_iwe)
									# 		for i in temp5:
									# 			delta = (i[1]-i[0]+1)
									# 			frame = delta/3
									# 			if delta%3==0:
									# 				print('delta=', delta, 'div', frame)
									# 			else:
									# 				print('delta=', delta, 'div', frame, q_gene_id, tracking_q_transcript_id, coordinates, q_strand)
									else:
										i_ic_longer_ref_subset_qry_refnotzero_other_c +=1				
								else:
									i_ic_longer_ref_subset_qry_other_c+=1					#    0
							else:
								i_ic_longer_ref_notsubset_qry_c+=1							#1,683
								start_pos_intron_chain = intron_chain[0]
								s_start_pos_intron_chain = intron_chain[0][0]
								end_pos_intron_chain = intron_chain[-1] 
								e_end_pos_intron_chain = intron_chain[-1][-1] 
								start_pos_ref_intron_chain = ref_intron_chain[0]
								s_start_pos_ref_intron_chain = ref_intron_chain[0][0]
								e_end_pos_ref_intron_chain= ref_intron_chain[-1][-1]
								if s_start_pos_intron_chain ==s_start_pos_ref_intron_chain and e_end_pos_intron_chain == e_end_pos_ref_intron_chain: #CHANGES ARE ALWAYS INTERNAL OF LEFTMOST AND RIGHTMOST OF THE TRANSCRIPT!
									if temp5[0][0]==temp6[0][0] and temp5[-1][-1]==temp6[0][-1]:
										if tracking_q_transcript_id not in isoform_dict.keys():
											isoform_dict[tracking_q_transcript_id]=1
											category = "Internal transcript boundaries, within intron" #not limited to 1 intron
											i_ic_longer_ref_notsubset_qry_internalwithinintron_c+=1
											#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
											#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = iwi_pid)
											output_features(q_transcript, o_iwi)
									else:
										if tracking_q_transcript_id not in isoform_dict.keys():
											isoform_dict[tracking_q_transcript_id]=1
											category = "Complex changes, internal transcript boundaries" #internal and/or other modifications
											i_ic_longer_ref_notsubset_qry_complexchanges_c+=1	# 167
											#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
											#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = complexchanges_pid)
											output_features(q_transcript, o_cc)
								else:
									build_linked_isoform(tracking_q_transcript_id, temp5)
									if tracking_q_transcript_id not in isoform_dict.keys():
										isoform_dict[tracking_q_transcript_id]=1
										if tracking_q_transcript_id in linked_dict_isoform.keys():
											category = "Linked"
											i_ic_longer_ref_subset_qry_other_linked_c+=1
											output_features(q_transcript, o_merged)
											k, v = tracking_q_transcript_id, linked_dict_isoform.get(tracking_q_transcript_id)
											#print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:], '2 genes', sep= '\t') if len(v)==5 else print(category, tracking_q_transcript_id, coordinates, q_strand, v[2:],'3 genes', sep= '\t')
											# > print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t') if len(v)==5 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t')
											#--print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:], '2 genes', sep= '\t', file = o_merged) if len(v)==4 else print(q_gene_id, tracking_q_transcript_id, coordinates, q_strand, category,v[2:],'3 genes', sep= '\t', file = o_merged)
										else:
											category = "Complex changes, diff transcript boundaries" #diff transcript boundaries
											i_ic_longer_ref_notsubset_qry_other_c+=1					#1,027
											#ori_PID, rev_PID, delta_PID = calc_rev_PID(q_transcript_id, coordinates, q_strand)
											#print(tracking_q_transcript_id, coordinates, q_strand, ori_PID, rev_PID, delta_PID, sep ='\t', file = complexchanges_pid)
											output_features(q_transcript, o_cc)
						elif len(intron_chain)<len(ref_intron_chain):
							#print('shorter', k, tracking_q_transcript_id, intron_chain, len(intron_chain), 'ref', ref_intron_chain, len(ref_intron_chain), ref_coordinate)
							if tracking_q_transcript_id not in isoform_dict.keys():
								isoform_dict[tracking_q_transcript_id]=1
								category = 'Shorter'
								i_ic_shorter_c+=1
								output_features(q_transcript, o_short)
		elif 'XLOC' in q_transcript_id:
			transcript_id = q_attributes['ID'][0]
			q_xloc_c+=1
			category = 'Novel'
			parent_novel(q_transcript, o_novel)
			#--parent_novel(q_transcript, o_overlap_novel)
			#--parent_novel(q_transcript, o_novel_only)
			#--print(q_transcript_id, transcript_id, coordinates, q_strand, category, sep= '\t', file = o_novel)

		#if intron_chain in ref_transcript_intron_dict.values() or intron_chain in ref_transcript_intron_isoform_dict.values():
		#	putative_novel_c+=1 #159
		#	category = "novel"
		#	for gene in genes:
		#		print (gene)
		#	print(q_transcript)
		#	for feature in features:
		#		print(feature)
		#	print(coordinates, q_strand, category, sep= '\t')#, file=n_out)

#print(c_conserved, file=sys.stderr)
# for k,v in sorted(linked_dict.items()):
# 	#print(k,v)
# 	if len(v)==4:
# 		print (k,v, '2 genes')
# 	else:
# 		print(k,v,'3 genes')
# print(len(linked_dict))

# for k,v in sorted(linked_dict_isoform.items()):
# 	#print(k,v)
# 	if len(v)==4:
# 		print (k,v, '2 genes')
# 	else:
# 		print(k,v,'3 genes')
# print(len(linked_dict_isoform))

print('', file=sys.stderr)

# print the results
print('Query path:\t{}'.format(qry))
print('Reference path:\t{}'.format(ref))
print('#'*120)
print('Reference:')
print('Total protein-coding genes:\t\t\t\t{}'.format(g_c))
print('Total protein-coding transcripts:\t\t\t{}'.format(t_c))
print('..transcripts with isoforms:\t\t\t {}'.format(isoform_c))
print('..transcripts with no isoforms:\t\t\t{}'.format(no_isoform_c))
print('#'*120)
print('Query:')
print('Total protein-coding transcripts:\t\t\t{}'.format(qt_c))
print('..non XLOC transcripts/transcripts with WBGene:\t\t{}'.format(q_notxloc_c))
print('....without isoform:\t\t\t\t{}'.format(q_noisoform_c))
print('......same length intron chain:\t\t{}'.format(ic_same_c))
print('      -no CDS transcript (excluded):\t{}\t{}'.format(ic_same_noCDStranscript_c, 'note: also exclude single exons with UTR intron'))
print('      -completely match:\t{}'.format(ic_same_completelymatch_c))
print('      -different content:\t{}'.format(ic_same_differentcontent_c))
print('      --single intron:\t{}'.format(ic_same_diff_singleintron_c))
print('      ---AA:\t{}'.format(ic_same_diff_singleintron_AA_c))
print('      ---AD:\t{}'.format(ic_same_diff_singleintron_AD_c))
print('      ---Both AA&AD:\t{}'.format(ic_same_diff_singleintron_both_AA_AD_c))
print('      --more than one intron:\t{}'.format(ic_same_diff_notsingleintron_c))
print('      ---1 difference:\t{}'.format(ic_same_diff_notsingleintron_1diff_c))
print('      ----AA:\t{}'.format(ic_same_diff_notsingleintron_1diff_AA_c))
print('      ----AD:\t{}'.format(ic_same_diff_notsingleintron_1diff_AD_c))
print('      ----Both AA&AD:\t{}'.format(ic_same_diff_notsingleintron_1diff_both_AA_AD_c))
print('      --->1 difference (same length, diff content, complex changes, multiple AS events):\t{}'.format(ic_same_diff_notsingleintron_morediff_c))
print('......longer intron chain:\t\t {}'.format(ic_longer_c))
print('      -ref is subset of qry:\t{}'.format(ic_longer_ref_subset_qry_c))
print('      --ref is subset, ref is zero (ref CDS absent, qry present):\t{}'.format(ic_longer_ref_subset_qry_refiszero_c))
print('      --ref is subset, ref not zero:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_c))
print('      ---5prime linked:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_left_linked_c))
print('      ---5prime extension:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_left_extension_c))
print('      ---3prime linked:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_right_linked_c))
print('      ---3prime extension:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_right_extension_c))
print('      ---5prime & 3prime new introns, linked:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_left_right_linked_c))
print('      ---both 5prime & 3prime extension:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c))
print('      ---complex changes (both 5 & 3prime extension+iwe):\t{}'.format(ic_longer_ref_subset_qry_refnotzero_complexchanges_ext_iwe_c))
print('      ---internal intron within exon:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c))
print('      ---other:\t{}'.format(ic_longer_ref_subset_qry_refnotzero_other_c))
print('      --ref is subset, other:\t{}'.format(ic_longer_ref_subset_qry_other_c))
print('      -ref is not subset of qry:\t{}'.format(ic_longer_ref_notsubset_qry_c))
print('      --ref is not subset, internal of transcript boundaries changes, within intron (>= 1 additional intron within intron):\t{}'.format(ic_longer_ref_notsubset_qry_internalwithinintron_c))
print('      --ref is not subset (complex changes, internal of transcript boundaries):\t{}'.format(ic_longer_ref_notsubset_qry_complexchanges_c))
print('      --ref is not subset, linked (merging and 1 or more diff junctions)\t{}'.format(ic_longer_ref_subset_qry_other_linked_c))
print('      --ref is not subset, other (complex changes, unrelated to transcript boundaries)\t{}'.format(ic_longer_ref_notsubset_qry_other_c))
print('......shorter intron chain (ignored):\t{}'.format(ic_shorter_c))

print('....with isoform:\t\t\t\t  {}'.format(len(isoform_dict)))

print('......same length intron chain:\t\t{}'.format(i_ic_same_noCDStranscript_c+i_ic_same_completelymatch_c+i_ic_same_diff_singleintron_AA_c+i_ic_same_diff_singleintron_AD_c+i_ic_same_diff_singleintron_both_AA_AD_c+i_ic_same_diff_notsingleintron_1diff_AA_c+i_ic_same_diff_notsingleintron_1diff_AD_c+i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c+i_ic_same_diff_notsingleintron_morediff_c)) #i_ic_same_c (contains dupl)
print('      -no CDS transcript (excluded):\t{}'.format(i_ic_same_noCDStranscript_c))
print('      -completely match:\t{}'.format(i_ic_same_completelymatch_c))
print('      -different content:\t{}'.format(i_ic_same_diff_singleintron_AA_c+i_ic_same_diff_singleintron_AD_c+i_ic_same_diff_singleintron_both_AA_AD_c+i_ic_same_diff_notsingleintron_1diff_AA_c+i_ic_same_diff_notsingleintron_1diff_AD_c+i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c+i_ic_same_diff_notsingleintron_morediff_c)) #i_ic_same_differentcontent_c (contains dupl)
print('      --single intron:\t{}'.format(i_ic_same_diff_singleintron_AA_c+i_ic_same_diff_singleintron_AD_c+i_ic_same_diff_singleintron_both_AA_AD_c)) #i_ic_same_diff_singleintron_c (contains dupl)
print('      ---AA:\t{}'.format(i_ic_same_diff_singleintron_AA_c))
print('      ---AD:\t{}'.format(i_ic_same_diff_singleintron_AD_c))
print('      ---Both AA&AD:\t{}'.format(i_ic_same_diff_singleintron_both_AA_AD_c))
print('      --more than one intron:\t{}'.format(i_ic_same_diff_notsingleintron_1diff_AA_c+i_ic_same_diff_notsingleintron_1diff_AD_c+i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c+i_ic_same_diff_notsingleintron_morediff_c)) #i_ic_same_diff_notsingleintron_c (contains dupl)
print('      ---1 difference:\t{}'.format(i_ic_same_diff_notsingleintron_1diff_AA_c+i_ic_same_diff_notsingleintron_1diff_AD_c+i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c)) #i_ic_same_diff_notsingleintron_1diff_c (contains dupl)
print('      ----AA: {}, AD: {}, both AA&AD: {}'.format(i_ic_same_diff_notsingleintron_1diff_AA_c, i_ic_same_diff_notsingleintron_1diff_AD_c, i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c))
print('      --->1 difference (same length, diff content, complex changes, multiple AS events):\t{}'.format(i_ic_same_diff_notsingleintron_morediff_c))

print('......longer intron chain:\t\t {}'.format(i_ic_longer_ref_subset_qry_refiszero_c+i_ic_longer_ref_subset_qry_refnotzero_left_extension_c+i_ic_longer_ref_subset_qry_refnotzero_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_left_linked_c+i_ic_longer_ref_subset_qry_refnotzero_right_linked_c+i_ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c+i_ic_longer_ref_subset_qry_refnotzero_other_c+i_ic_longer_ref_subset_qry_other_c+i_ic_longer_ref_notsubset_qry_internalwithinintron_c+i_ic_longer_ref_notsubset_qry_complexchanges_c+i_ic_longer_ref_subset_qry_other_linked_c+i_ic_longer_ref_notsubset_qry_other_c)) #i_ic_longer_c(contains dupl)
print('      -ref is subset of qry:\t{}'.format(i_ic_longer_ref_subset_qry_refiszero_c+i_ic_longer_ref_subset_qry_refnotzero_left_extension_c+i_ic_longer_ref_subset_qry_refnotzero_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_left_linked_c+i_ic_longer_ref_subset_qry_refnotzero_right_linked_c+i_ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c+i_ic_longer_ref_subset_qry_refnotzero_other_c+i_ic_longer_ref_subset_qry_other_c)) #i_ic_longer_ref_subset_qry_c (contains dupl)
print('      --ref is subset, ref is zero:\t{}'.format(i_ic_longer_ref_subset_qry_refiszero_c))
print('      --ref is subset, ref not zero:\t{}'.format(i_ic_longer_ref_subset_qry_refnotzero_left_extension_c+i_ic_longer_ref_subset_qry_refnotzero_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_left_linked_c+i_ic_longer_ref_subset_qry_refnotzero_right_linked_c+i_ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c+i_ic_longer_ref_subset_qry_refnotzero_other_c)) #i_ic_longer_ref_subset_qry_refnotzero_c (contains dupl)
print('      ---5prime extension: {}, 3prime extension: {}, both 5prime & 3prime extension: {}'.format(i_ic_longer_ref_subset_qry_refnotzero_left_extension_c, i_ic_longer_ref_subset_qry_refnotzero_right_extension_c, i_ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c))
print('      ---5prime linked: {}, 3prime linked: {}, both 5prime & 3prime linked: {}'.format(i_ic_longer_ref_subset_qry_refnotzero_left_linked_c, i_ic_longer_ref_subset_qry_refnotzero_right_linked_c, 0))
print('      ---internal intron within exon:\t{}'.format(i_ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c))
print('      ---other:\t{}'.format(i_ic_longer_ref_subset_qry_refnotzero_other_c))
print('      --ref is subset, other:\t{}'.format(i_ic_longer_ref_subset_qry_other_c))
print('      -ref is not subset of qry:\t{}'.format(i_ic_longer_ref_notsubset_qry_internalwithinintron_c+i_ic_longer_ref_notsubset_qry_complexchanges_c+i_ic_longer_ref_subset_qry_other_linked_c+i_ic_longer_ref_notsubset_qry_other_c)) #i_ic_longer_ref_notsubset_qry_c (contains dupl)
print('      --ref is not subset, internal of transcript boundaries changes, within intron (>= 1 additional intron within intron):\t{}'.format(i_ic_longer_ref_notsubset_qry_internalwithinintron_c))
print('      --ref is not subset (complex changes, internal of transcript boundaries):\t{}'.format(i_ic_longer_ref_notsubset_qry_complexchanges_c))
print('      --ref is not subset, linked (merging and 1 or more diff junctions)\t{}'.format(i_ic_longer_ref_subset_qry_other_linked_c))
print('      --ref is not subset, other (complex changes, unrelated to transcript boundaries)\t{}'.format(i_ic_longer_ref_notsubset_qry_other_c))

print('......shorter intron chain (ignored), note some are possibly useful:\t{}'.format(i_ic_shorter_c))

print('..XLOC transcripts/novel:\t\t\t\t  {}'.format(q_xloc_c))
print('\n')
print('='*150)
print('Recap based on categories:')
print("Category 1. Complete match/WB confirmed \t\t{}".format(ic_same_completelymatch_c+i_ic_same_completelymatch_c))
print("Category 2. 3' extension \t\t\t\t {}".format(ic_longer_ref_subset_qry_refnotzero_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_right_extension_c))
print("Category 3. 5' extension \t\t\t\t {}".format(ic_longer_ref_subset_qry_refnotzero_left_extension_c+i_ic_longer_ref_subset_qry_refnotzero_left_extension_c))
print("Category 4. 5' & 3' extension \t\t\t\t  {}".format(ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c+i_ic_longer_ref_subset_qry_refnotzero_both_left_right_extension_c))
print("Category 5. Contain intron overlapping internal exon \t {}".format(ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c+i_ic_longer_ref_subset_qry_refnotzero_internalwithinexon_c))
print("Category 6. Contain intron overlapping intron \t\t {}".format(ic_longer_ref_notsubset_qry_internalwithinintron_c+i_ic_longer_ref_notsubset_qry_internalwithinintron_c))
print("Category 7. Alternative donor \t\t\t\t {}".format(ic_same_diff_singleintron_AD_c+ic_same_diff_notsingleintron_1diff_AD_c+i_ic_same_diff_singleintron_AD_c+i_ic_same_diff_notsingleintron_1diff_AD_c))
print("Category 8. Alternative acceptor \t\t\t {}".format(ic_same_diff_singleintron_AA_c+ic_same_diff_notsingleintron_1diff_AA_c+i_ic_same_diff_singleintron_AA_c+i_ic_same_diff_notsingleintron_1diff_AA_c))
print("Category 9. Alternative donor & acceptor \t\t {}".format(ic_same_diff_singleintron_both_AA_AD_c+ic_same_diff_notsingleintron_1diff_both_AA_AD_c+i_ic_same_diff_singleintron_both_AA_AD_c+i_ic_same_diff_notsingleintron_1diff_both_AA_AD_c))
print("Category 10. Merging 2 or more genes \t\t\t {}".format(ic_longer_ref_subset_qry_refnotzero_left_linked_c+ic_longer_ref_subset_qry_refnotzero_right_linked_c+ic_longer_ref_subset_qry_refnotzero_left_right_linked_c+ic_longer_ref_subset_qry_other_linked_c+i_ic_longer_ref_subset_qry_refnotzero_left_linked_c+i_ic_longer_ref_subset_qry_refnotzero_right_linked_c+0))
print("Category 11. Complex changes \t\t\t\t{}".format(ic_same_diff_notsingleintron_morediff_c+ic_longer_ref_subset_qry_refnotzero_complexchanges_ext_iwe_c+ic_longer_ref_notsubset_qry_complexchanges_c+ic_longer_ref_notsubset_qry_other_c+i_ic_same_diff_notsingleintron_morediff_c+i_ic_longer_ref_notsubset_qry_complexchanges_c+i_ic_longer_ref_notsubset_qry_other_c))
print("Category 12. Novel \t\t\t\t\t {}".format(q_xloc_c))
print('\n')
print('Other categories')
print("Category 13. Single-intron transcripts in both/only in ref {}".format(ic_same_noCDStranscript_c+i_ic_same_noCDStranscript_c+ic_longer_ref_subset_qry_refiszero_c+i_ic_longer_ref_subset_qry_refiszero_c))
print("Category 14. Shorter \t\t\t\t\t{}".format(ic_shorter_c+i_ic_shorter_c))
print('\n')
print('Done!', file=sys.stderr)
