#import sys, pysam, pprint, argparse
import sys, argparse, pprint, os, errno, requests, json, re, copy, tempfile, pysam
from statistics import mean
#import sys, argparse, pprint, os, errno, requests, json, re, copy, tempfile
from pyfaidx import Faidx
from Bio.Seq import Seq

def get_all_mapping_qualities_for_one_locus(samfile, query_chr, chr_start, chr_stop):
	chr_start = int(chr_start)
	chr_stop  = int(chr_stop)
	#print("%s,%s,%s")%(query_chr,chr_start,chr_stop)
	mapq_list = []
	for read in samfile.fetch(query_chr, chr_start - 1, chr_stop):
		#print(read.mapq)
		mapq_list.append(read.mapq)
	return(mean(mapq_list))

def tally_non_ref_one_position(query_chr, query_chr_pos, samfile):
	ref_base_from_fasta = FA.fetch(query_chr, query_chr_pos, query_chr_pos).seq
	non_ref_total = 0
	total_cov     = 0
	##NOTE: PCR duplicates are silently filtered out by pileup.  if reads are lower counts than expected this may be why.
	for pileupcolumn in samfile.pileup(query_chr, query_chr_pos - 1, query_chr_pos):
		if pileupcolumn.pos == query_chr_pos - 1:
			#print ("\ncoverage at base %s,%s = %s" %(pileupcolumn.pos + 1, query_chr, pileupcolumn.n))
			for pileupread in pileupcolumn.pileups:
				if pileupread.query_position is None:
					non_ref_total += 1
				else:
					read_call = pileupread.alignment.query_sequence[pileupread.query_position]
					if ref_base_from_fasta.upper() != read_call.upper():
						non_ref_total += 1
				total_cov += 1
	ref_total = total_cov - non_ref_total
	#print(("%s,%s,%s")%(non_ref_total,total_cov,ref_total))
	return total_cov, non_ref_total, ref_total

def evaluate_criteria_5_for_one_locus(chr,chr_pos,flank_start,flank_stop,samfile,exclude_array):
	fail = False
	total_positions_gte_percent_cutoff = 0
	failure_pos = []
	for x in range(flank_start, flank_stop):
		total_cov, non_ref_total, ref_total = tally_non_ref_one_position(query_chr=chr,query_chr_pos=x,samfile=samfile)
		non_ref_freq = 0.0
		if total_cov > 0:
			non_ref_freq = non_ref_total/float(total_cov)
		if total_cov >= ARGS['min_percent_cov'] and non_ref_freq > ARGS['non_ref_freq'] and x not in exclude_array:
			total_positions_gte_percent_cutoff += 1
			failure_pos.append(x)
	if total_positions_gte_percent_cutoff >= ARGS['pos_total']:
		fail = True
	return fail,total_positions_gte_percent_cutoff,failure_pos

def process_vcf_file(vcf_file, bam_file):
	samfile = pysam.AlignmentFile(bam_file, "rb")
	for line in open(vcf_file):
		line = line.rstrip()
		if line[0] == '#':
			continue
		line = line.split("\t")
		chr, chr_pos = line[:2]
		ref_allele, var_allele = line[3:5]
		ref_base = ref_allele[0]
		chr_pos = int(chr_pos)
		flank_start = chr_pos - ARGS['flanking_bp']
		flank_stop  = chr_pos + ARGS['flanking_bp']

		variant_max_allele_length = len(ref_allele)
		if len(var_allele) > variant_max_allele_length:
			variant_max_allele_length = len(var_allele)

		var_id = chr + '@' + str(chr_pos) + '|' + ref_allele + '|' + var_allele

		exclude_start = chr_pos - variant_max_allele_length + 1
		exclude_stop  = chr_pos + variant_max_allele_length - 1
		exclude_array = []
		for i in range(exclude_start, exclude_stop + 1):
			exclude_array.append(i)

		exclude_array = []

		#print((var_id + "\t%s\t%s,%s,%s")%(variant_max_allele_length,exclude_start,exclude_stop,exclude_array))hr
		#continue

		#ref_fasta = FA.fetch(chr, chr_pos, chr_pos).seq #sanity check- looks good (agrees w/vcf)
		#print("%s,%s,%s,%s,%s")%(chr,chr_pos,chr_pos,ref_base, ref_fasta)
		fail,total_positions_gte_percent_cutoff,failure_pos = evaluate_criteria_5_for_one_locus(chr=chr, chr_pos=chr_pos, flank_start=flank_start, flank_stop=flank_stop, samfile=samfile, exclude_array=exclude_array)
		mapping_qual_mean = get_all_mapping_qualities_for_one_locus(samfile=samfile,query_chr=chr,chr_start=chr_pos,chr_stop=chr_pos)
		print((var_id + "\t" + "%s,%s,%s,avg_mapping_qual=%s")%(fail,total_positions_gte_percent_cutoff,failure_pos,mapping_qual_mean))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--flanking_bp',     default=70,   help='+/- bp window surrounding the candidate indel.', type=int)
	parser.add_argument('-p', '--pos_total',       default=5,    help="positions that have non-reference base frequencies >Xpercent.", type=int)
	parser.add_argument('-n', '--non_ref_freq',    default=0.02, help='non-reference base frequency to use as threshold.  i.e. 0.02=2percent', type=float)
	parser.add_argument('-m', '--min_percent_cov', default=20,   help='the minimum coverage needed to consider the non_ref_freq for a position.', type=int)
	parser.add_argument('-v', '--vcf_file', help='VCF file to operate on.', type=str)
	parser.add_argument('-b', '--bam_file', help='BAM file to operate on.', type=str)
	parser.add_argument('-r', '--ref',      help='indexed reference in FASTA format i.e. human_g1k_v37_decoy.fasta', type=str)
	parser.add_argument('--query_chr',      help='query this chromosome (ignore VCF file)', type=str)
	parser.add_argument('--query_chr_pos',  help='query this chromosome position (ignore VCF file)', type=int)
	args = parser.parse_args()
	args = vars(args)
	ARGS = args
	if len(sys.argv) == 1: 
		parser.print_help()
		sys.exit() 
	FA = Faidx(ARGS['ref'])
	if 'vcf_file' in ARGS:
		process_vcf_file(args['vcf_file'], args['bam_file'])
	elif 'query_chr' in ARGS and query_chr_pos in ARGS:
		print("Not implemented yet.")

