#!/usr/bin/env/python
# -*- coding: utf-8 -*- 

__author__ = 'xxx'
__version__ = '0.1.0'
__date__ = 'Sep 2017'

"""
Strain-GeMS: identifying strains with metagenomic reads

Copyright(c) 2017 xxx)

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

https:xxx



USAGE = \
"""Usage: 
%prog [options] -c/--conf <config.file> -o/--outdir <output directory>

The config file should follow the format:

//
sample: [sample1_ID]
fq1: [forward reads fastq]
fq2: [reverse/mate reads fastq]
metaphlan2: [metaphlan2 output]  
//
sample: [sample2_ID]
fq1: [forward reads fastq]
fq2: [reverse/mate reads fastq]
metaphlan2: [metaphlan2 output]  
...


Or, if only one fastq file per sample, the config file should look like this:

//
sample: [sample1_ID]
fq: [sample reads fastq]
//
sample: [sample2_ID]
fq: [sample reads fastq]
...

Add --help to see a full list of required arguments to run Strain-GeMS.



Strain-GeMS: identifying strains with metagenomic reads

Strain-GeMS aims to profile the strains within a species using short read metagenomes.
It deploys SNP patterns based on statistically sound SNP calling and optimzied SNP clustering 
to recover strain genotypes within a species using a set of genes.


Copyright: xxx
"""

import sys, os, re, glob, random, shutil, resource
from copy import copy
from optparse import OptionParser, OptionGroup
from itertools import combinations, product, permutations
from operator import itemgetter, attrgetter
import multiprocessing as mp
from subprocess import PIPE, Popen, call
import tarfile, bz2, cPickle

import numpy as np
from numpy.random import uniform, binomial, dirichlet, normal
from scipy.spatial.distance import cdist, pdist, euclidean
from scipy.cluster import hierarchy, vq
from scipy.stats import percentileofscore, entropy
import networkx as nx
from Bio import SeqIO

MAX_ITER = 10000
TOL = 0.01
SIGMA = 0.1

class Sample:
	def __init__(self):
		self.name = None
		self.fq1 = None
		self.fq2 = None
		self.fq = None
		self.metaphlan2 = None
	
	def is_empty(self):
		if self.fq1 == None and self.fq2 == None and self.fq == None:
			return True
		if self.fq != None: return False
		if self.fq1 != None and self.fq2 != None: return False
		
	def get_files(self):
		fqs = []
		if self.fq != None: fqs.append(self.fq)
		elif self.fq1 != None and self.fq2 != None: fqs = [self.fq1, self.fq2]
		return fqs, self.metaphlan2
		
	def clear(self):
		self.name = None
		self.fq1 = None
		self.fq2 = None
		self.fq = None
		self.metaphlan2 = None
		
class Pileup:
	def __init__(self):
		self.pid = None
		self.pos = 0
		self.ref_base = None
		self.cov = []
		self.pileup = []
		
	def update(self, index, cov, p):
		self.cov[index] = cov
		self.pileup[index] = p

class SNP:
	def __init__(self):
		self.gene = None
		self.loc = 0
		self.counts = []
		self.freqs = []
		self.barcode = []
		self.sorted_freqs = []
		self.sorted_counts = []
		self.sorted_bases = []
	
	def counts2freqs(self, minor_cov):
		freqs = []
		for count in self.counts:
			f = [-1, -1, -1, -1]
			if sum(count) < minor_cov: pass
			else:
				for i in range(4):
					f[i] = float(count[i])/sum(count)
			freqs.append(f)
		self.freqs = freqs
	
	def sort_freqs(self):
		for count, freq in zip(self.counts, self.freqs):
			x = []
			for f, c, b in zip(freq, count, ['A','C','G','T']): x.append((f, c, b))
			self.sorted_freqs.append(map(itemgetter(0), sorted(x, key = lambda a: a[0], reverse = True)))
			self.sorted_counts.append(map(itemgetter(1), sorted(x, key = lambda a: a[0], reverse = True)))
			self.sorted_bases.append(map(itemgetter(2), sorted(x, key = lambda a: a[0], reverse = True)))
			
	def total_count(self):
		x = 0
		for counts in self.counts: x+=sum(counts)
		return x
	
	def loci(self, minor_cov):
		loci = []
		N = len(self.counts)
		for sampleIndex in range(N):
			for b, c in zip(self.sorted_bases[sampleIndex], self.sorted_counts[sampleIndex]):
				if c >= minor_cov and b not in loci: loci.append(b)
		return loci

class Model:
	def __init__(self):
		self.profile = []
		self.err = []
		self.overall_err = 0
		self.numSNPs = 0
		self.locs = []
		self.strains = []
		self.aligned_locs = []
		self.aligned_strains = []
		self.AICc = 1e5
		self.RSS = 1e5
		self.AIC = 1e5
		
	def numStr(self):
		return len(self.strains)
	
	def effective_num_SNPs(self):
		return len(self.strains[0]) - self.strains[0].count('-')
	
	def evaluate_model(self, SNPs, sample_index):
		if self.err == []:
			N = len(SNPs[0].counts)
			for i in range(N): self.err.append(0.0)
		counts = []
		dict_SNPs = {}
		for S in SNPs: dict_SNPs[(S.gene, S.loc)] = S
		for i in range(len(self.locs)):
			S = dict_SNPs[tuple(self.locs[i])]
			counts.append(S.counts[sample_index])
		M, O, F = to_matrices(self.strains, counts)
		self.err[sample_index] = fit_err(np.array(self.profile), M, F)	
	
	def align_strains(self, locs):
		self.aligned_locs = locs
		dict_SNPs = {}
		for loc_index in range(len(self.locs)):
			dict_SNPs[tuple(self.locs[loc_index])] = loc_index
		L = len(self.locs)
		N = len(self.strains)
		for i in range(N): self.aligned_strains.append('')
		for loc in locs:
			if loc not in dict_SNPs:
				for i in range(N): self.aligned_strains[i] += '-'
				continue
			else:
				loc_index = dict_SNPs[loc]
				for i in range(N): self.aligned_strains[i] += self.strains[i][loc_index]
			
	def calAICc(self, numStr, sample_index): # using the alternative AIC by RSS
		N = numStr-1
		numSNPs = self.effective_num_SNPs()
		corr_factor = float(2*N*(N+1))/(numSNPs-N-1)
		RSS = self.err[sample_index] * numSNPs
		AIC = N*np.log(RSS/N) + 2*N
		self.AIC = AIC
		self.RSS = RSS
		self.AICc = AIC + corr_factor

class SNP_layer:
	def __init__(self, SNPs, genotypes, profiles, w, err):
		self.genotypes = genotypes
		self.profiles = profiles
		self.SNPs = SNPs
		self.weight = w
		self.err = err
		
############################ FUNCTIONS ##############################
def index_base(index):
	if index == 0: return 'A'
	if index == 1: return 'C'
	if index == 2: return 'G'
	if index == 3: return 'T'

def base_index(base):
	if base == 'A' or base == 'a': return 0
	if base == 'C' or base == 'c': return 1
	if base == 'G' or base == 'g': return 2
	if base == 'T' or base == 't': return 3
	return 0	
# End of base_index

def is_snp(counts):
	max_depth = 0
	flag = 0
	for count in counts:
		if sum(count) > max_depth: max_depth = sum(count)
		sorted_count = sorted(count, reverse = True)
		if sum(sorted_count[2:]) > 0: return False
		if sum(count) < 10: continue
		if sorted_count[1] > 2 and float(sorted_count[1])/sum(count) > 0.05:
			flag = 1
	if max_depth < 10: flag = 0
	if flag == 0: return False
	return True

def freq_err(x, y, guided = False):
	err = 0.
	if guided:
		for a, b in zip(x, y):
			err += abs(a-b)
	else:
		# knapsack problem using dynamic programming
		# find sums of subsets of y so err(x, y) is minimized
		# now just a brute-force solution for fast implementation
		sorted_x = sorted(x, reverse = True)[:2]
		splits = []
		for l in range(1, len(y)):
			for sub in list(combinations(y, l)):
				substract = [a for a in y if a not in sub]
				combo = [list(sub), substract]
				splits.append(combo)
		min_err = 1
		for split in splits:
			f = []
			e = 0.
			for chunk in split:
				f.append(sum(chunk))
			f = sorted(f, reverse = True)
			for a, b in zip(sorted_x, f):
				e += abs(a-b)
			if e < min_err: min_err = e
		err = min_err
	return err


def decode_base(ref_base, pstring):
	res = interpret_pstring(ref_base, pstring)
	base = sorted(zip(['A','C','G','T'], res['bases']), key = lambda x: x[1], reverse = True)[0][0]
	return base


def pileup2counts(pileup, ref_base):
	x = [0,0,0,0]
	ref_index = base_index(ref_base)
	x[ref_index] += pileup.count(',')+pileup.count('.')
	for i in range(4):
		x[i] += pileup.count(index_base(i)) + pileup.count(index_base(i).lower())
	return x
	

def construct_snpflow(options, samples, species_list):
	samtools = options.samtools
	N = len(samples)
	num_proc = options.num_proc
	quiet = options.quiet
	nofile_limit = resource.RLIMIT_NOFILE
	pileup_dir = options.outdir+'/pileups/'
	sam_dir = options.outdir + '/sams/'
	pileup_dir = options.outdir + '/pileups/'
	ref_file = sam_dir + 'merged_ref.ffn'
	if not os.path.exists(ref_file):
		ref_file = options.outdir+'/merged_ref.ffn'
	
	ref_seqs = {}
	for record in SeqIO.parse(ref_file, 'fasta'):
		species = record.name.split('|')[-1].replace('s__','')
		pid = record.name.split('|')[-2]
		if species not in ref_seqs: ref_seqs[species] = {}
		ref_seqs[species][pid] = (record.name, record.seq)
	
	if len(glob.glob(pileup_dir+'*.pileups')) == 0 or (not os.path.exists(pileup_dir)):
		if not quiet: sys.stdout.write('  Now mpileup reads...\n')
		bamfiles = []
		mpileup_cmds = []
		for sample in samples:
			bamfile = sam_dir + 'reads_vs_ref.'+sample.name+'.bam'
			bamfiles.append(os.path.realpath(bamfile))
		# mpileup of all bams
		mpileup_file = sam_dir + 'reads_vs_ref.mpileup'
		if not os.path.exists(mpileup_file):
			mfh = open(mpileup_file, 'w')
			logfile = mpileup_file+'.log'
			logfh = open(logfile, 'w')
			#mpileup_cmd = [samtools, 'mpileup', '-C50', '-d10000', '-Df', ref_file]
			mpileup_cmd = [samtools, 'mpileup','-s','-C50', '-d10000', '-Df', ref_file]
			for bamfile in bamfiles: mpileup_cmd.append(bamfile)
			p1=Popen(mpileup_cmd, stdout = mfh, stderr = logfh)
			p1_stdout, p1_stderr = p1.communicate()
			mfh.close()
			logfh.close()
		if not quiet: sys.stdout.write('  Done. Now splitting into per species mpileups...\n')
		
		# split the mpileup file into species files
		for species in ref_seqs:
			sp_split_logfh = open(sam_dir+species+'.mpileup.log', 'w')
			vcf_file = sam_dir + species+'.mpileup'
			vfh = open(vcf_file, 'w')
			p1 = Popen(['cat', mpileup_file], stdout = PIPE, stderr = sp_split_logfh)
			p2 = Popen(['grep', species], stdin = p1.stdout, stdout = vfh, stderr = sp_split_logfh)
			p1.stdout.close()
			output, err = p2.communicate()
			sp_split_logfh.close()
			vfh.close()
			if not quiet: sys.stdout.write('    %s\n' % species)
		if not quiet: sys.stdout.write('  Done. Now formatting them...\n')
		
		# format mpileup file	
		if not os.path.exists(pileup_dir): os.mkdir(pileup_dir)
		for species in ref_seqs:
			vcf_file = sam_dir + species+'.mpileup'
			if not os.path.exists(vcf_file): continue
			if os.stat(vcf_file).st_size == 0:
				os.unlink(vcf_file)
				continue
			outfile = pileup_dir + '/' + species + '.pileups'
			logfile = outfile + '.log'
			gems_cmd = ['multigems', '-S', str(N), '-i', vcf_file, '-o', outfile]
			lfh = open(logfile, 'w')
			p1 = Popen(gems_cmd, stdout = lfh, stderr = lfh)
			p1_stdout, p1_stderr = p1.communicate()
			lfh.close()
			#ofh = open(outfile, 'w')
			#for line in open(vcf_file, 'r'):
				#cols = line[:-1].split('\t')
				#gene_id, pos, ref_base = cols[:3]
				#pid = gene_id.split('|')[-2]
				#covs = []
				#pileups = []
				#for x_index in range(3, len(cols), 3):
					#cov_st, pileup_st, q_st = cols[x_index:x_index+3]
					#covs.append(cov_st)
					#pileups.append(pileup_st)
				#ofh.write('%s\t%s\t%s\t%s\t%s\n' % (pid, pos, ref_base, ','.join(covs), '|'.join(pileups)))
			#ofh.close()
	if not quiet: sys.stdout.write('  Done!\n')
	
	# load ref sequences and call SNPs
	if not quiet: sys.stdout.write('  Now calling SNPs.\n')
	# locate all SNPs
	snpflow_dir = options.outdir + '/snpflows/'
	if not os.path.exists(snpflow_dir): os.mkdir(snpflow_dir)
	locateSNPs_cmds = []
	for sp in sorted(species_list):
		if sp.count('unclassified') == 1:
			if not quiet: sys.stdout.write('    [%s] Skipping because it is unclassified and potentially mixture of species.\n' % sp)
			continue
		pileup_file = pileup_dir+sp+'.pileups'
		snp_file = snpflow_dir+sp+'.snps'
		if not os.path.exists(pileup_file) or os.stat(pileup_file).st_size == 0:
			if not quiet: sys.stdout.write('    [%s] Skipping due to inadequate pileups.\n' % sp)
			continue
		if os.path.exists(snp_file) and os.stat(snp_file).st_size > 0:
			if not quiet: sys.stdout.write('    [%s] Previously finished.\n' % sp)
			continue
		locateSNPs_cmds.append([sp, pileup_file, ref_seqs[sp], snp_file, options])
	
	if not quiet:
		sp_to_process = map(itemgetter(0), locateSNPs_cmds)
		if not quiet:
			if len(sp_to_process) == 0:
				sys.stdout.write('  Done.\n')
				return 0
			sys.stdout.write('  %i species yet to process for SNPs.\n' % len(sp_to_process))
			for sp in sp_to_process: sys.stdout.write('    +%s\n' % sp)
		
	# mp mode
	if len(locateSNPs_cmds) > 0:
		pool = mp.Pool(options.num_proc)
		pool.map_async(locateSNPs, locateSNPs_cmds)
		pool.close()
		pool.join()
	
# End of construct_snpflow


################################ FUNCTIONS ####################################
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def need_metaphlan2(samples, outdir):
	need_mpa = False
	for sample in samples:
		if sample.metaphlan2 == None or not os.path.exists(sample.metaphlan2):
			need_mpa = True
			break
	return need_mpa
	"""
	if not os.path.exists(outdir): return need_mpa
	elif need_mpa == True:
		metaphlan2_dir = outdir + '/metaphlan2/'
		tmp_mpa = False
		for sample in samples:
			sample_metaphlan2 = metaphlan2_dir+sample.name+'.metaphlan2_out'
			if not os.path.exists(sample_metaphlan2):
				tmp_mpa = True
				break
		return tmp_mpa
	"""
	
# End of which

def run_metaphlan2(x):
	metaphlan2_py, mpa_pkl, bowtie2db, temp_dir, bowtie2out, metaphlan2_out, infile, log = x

	cmd = ["python", metaphlan2_py, "--mpa_pkl", mpa_pkl, "--input_type", "multifastq", "--nproc", "1", 
			"--tax_lev", "a", "--bowtie2db", bowtie2db, "-t", "rel_ab",
			"--bowtie2out", bowtie2out, "--bt2_ps", "very-sensitive",
			"--tmp_dir", temp_dir, "-o", metaphlan2_out, infile]
	
	logfh = open(log, 'w')
	p1 = Popen(cmd,stdout = logfh, stderr = logfh)
	p1_stdout, p1_stderr = p1.communicate()
#	call(cmd, stdout =sys.stdout, stderr = sys.stderr)
	logfh.close()
	
# End of run_metaphlan

def make_sp_list(samples, options):
	gsize_db = options.gsize_db
	min_cov = options.min_cov
	
	species_ab = {}
	species_cov = {}
	num_samples = len(samples)
	
	gsize = {}
	for line in open(gsize_db, 'rb'):
		clade, s_st, type_st = line[:-1].split('\t')
		s = float(s_st)
		if type_st == 'NA': continue
		gsize[clade] = float(s_st)
	
	for i, sample in enumerate(samples):
		mtp_out = sample.metaphlan2
		if sample.fq1 != None: read_file = sample.fq1
		else: read_file = sample.fq
		total_bases = 0
		extension = read_file.split('.')[-1]
		x_ext = 'fq'
		if extension == 'bz2': 
			x_ext = read_file.split('.')[-2]
			p1 = Popen(['bzcat', read_file], stdout = PIPE, stderr = sys.stderr)
		else:
			x_ext = read_file.split('.')[-1]
			p1 = Popen(['cat', read_file], stdout = PIPE, stderr = sys.stderr)
		p2 = Popen(['wc'], stdin = p1.stdout, stdout = PIPE, stderr = sys.stderr)
		p1.stdout.close()
		pipeout, pipeerr = p2.communicate()
		num_lines = int([a for a in pipeout.split(' ') if a != ''][0])
		
		# peek read length
		rlen = 0
		if extension == 'bz2':
			with bz2.BZ2File(read_file, 'rb') as ifh:
				tag = ifh.readline()
				seq = ifh.readline().rstrip('\n')
				rlen = len(seq)
		else:
			with open(read_file, 'rb') as ifh:
				tag = ifh.readline()
				seq = ifh.readline().rstrip('\n')
				rlen = len(seq)
		if x_ext == 'fa': total_bases = num_lines*rlen/2
		else: total_bases = num_lines*rlen/4
		if sample.fq1 != None: total_bases *= 2	
		
		for line in open(mtp_out, 'r'):
			if line[0] == '#': continue
			cols = line[:-1].split('\t')
			taxo = cols[0].split('|')
			if taxo[-1][:3] != 's__': continue
			if taxo[0][3:] == 'Viruses': continue
			sp_full = cols[0]
			sp = taxo[-1][3:]
			if sp.count('unclassified') > 0: continue # skip unclassified species in mpl2
			sp_gsize = None
			ab = float(cols[1])
			
			try: sp_gsize = gsize[sp_full]
			except KeyError:
				sys.stderr.write('[Warning]: Cannot locate the avg. genome size for species: %s, assuming it is ~4Mbp.\n' % sp)
				sp_gsize = 4
			
			cov = ab*(float(total_bases)/(sp_gsize*1e8))
			if sp not in species_ab: species_ab[sp] = [0.]*num_samples
			species_ab[sp][i] = ab
			if sp not in species_cov: species_cov[sp] = [0.]*num_samples
			species_cov[sp][i] = cov
			
	return species_ab, species_cov

# End of read_metaphlan_files

def metaphlan2_filter(samples, options):
	metaphlan2 = options.metaphlan2
	num_proc = options.num_proc
	outdir = options.outdir
	
	sample_names = []
	infiles = []
	metaphlan2_files = []
	num_samples = len(samples)
	for sample in samples:
		metaphlan2_files.append(sample.metaphlan2)
		if sample.fq1 != None: infiles.append((sample.fq1, sample.fq2))
		else: infiles.append((sample.fq))
	num_nones = len([x for x in metaphlan2_files if x == None])
	
	sp_list = {}
	if num_nones > 0:
		metaphlan2_dir = outdir + '/metaphlan2'
		if not os.path.exists(metaphlan2_dir): os.mkdir(metaphlan2_dir)
		temp_dir = metaphlan2_dir + '/tmp/'
		if not os.path.exists(temp_dir): os.mkdir(temp_dir)
		metaphlan2_py = os.path.abspath(which(metaphlan2))
		bowtie2dbs = glob.glob(os.path.dirname(metaphlan2_py)+'/db_v*/mpa_v*.bt2')
		if len(bowtie2dbs) == 0:
			sys.stderr.write('[FATAL]: Failure in locating bowtie2 db for MetaPhlAn2.\n')
			exit(1)
		bowtie2db = os.path.abspath(bowtie2dbs[0]).split('.')[0]
		try: mpa_pkl = glob.glob(os.path.dirname(metaphlan2_py)+'/db_v*/*.pkl')[0]
		except:
			sys.stderr.write('[FATAL]: Failure in locating mpa.pkl for MetaPhlAn2.\n')
			exit(1)
	
		# generate metaphlan2 profiles for those needed
		bowtie2_out = []
		metaphlan2_out = []
		input_file = []
		log_file = []
		for sample in samples:
			metaphlan2_outfile = metaphlan2_dir + '/' + sample.name + '.metaphlan2_out'
			if os.path.exists(metaphlan2_outfile): continue
			if sample.metaphlan2 != None: 
				p1 = Popen(["cp", sample.metaphlan2, metaphlan2_outfile])
				p1_stdout, p1_stderr = p1.communicate()
			else:
				metaphlan2_out.append(metaphlan2_outfile)
				bowtie2_out.append(temp_dir + '/bowtie2out.' + sample.name)
				log_file.append(temp_dir + '/metaphlan_log.' + sample.name)
				if sample.fq1 != None: infile = sample.fq1
				else: infile = sample.fq
				input_file.append(infile)
	
		if not options.quiet:
			sys.stdout.write("Profiling metagenome with MetaPhlAn2...\n")
		cmds = [[metaphlan2_py, mpa_pkl, bowtie2db, temp_dir, bowtie2out, metaphlan2_res, infile, log]
			for bowtie2out, metaphlan2_res, infile, log in 
			zip(bowtie2_out, metaphlan2_out, input_file, log_file)]
		
		#for cmd in cmds: run_metaphlan2(cmd)
		pool = mp.Pool(num_proc)
		pool.map_async(run_metaphlan2, cmds)
		pool.close()
		pool.join()
		# update sample.metaphlan2
		for sample in samples:
			metaphlan2_outfile = metaphlan2_dir + '/' + sample.name + '.metaphlan2_out'
			sample.metaphlan2 = metaphlan2_outfile
		if not options.quiet:
			sys.stdout.write("MetaPhlAn2 profiling done.\n")
		
		# check the run
		for log in log_file:
			for line in open(log, 'r'):
				if line.count('Error') > 0:
					sys.stderr.write('Error in running MetaPhlAn2:\n%s\n' % line)
					exit(1)
		shutil.rmtree(temp_dir)
	
	# select species:
	species_ab, species_cov = make_sp_list(samples, options)
	
	return species_ab, species_cov

# End of metaphlan_filter

def extract_tarball(db_tarfile, db_tarindex, bz2_dir, files_to_extract):
	with open(db_tarfile, 'rb') as tar:
		for line in open(db_tarindex, 'r'):
			m = line[:-1].rsplit(" ", 2)
			if m[0] in files_to_extract:
				tar.seek(int(m[1]))
				bz2chunk = tar.read(int(m[2]))
				bz2file = bz2_dir + m[0].replace('ref_db', '')
				bz2fh = open(bz2file, 'w')
				bz2fh.write(bz2chunk)
				bz2fh.close()
				
# End of extract_tarball

def make_refDB(species_list, db_tarfile, db_tarindex, options):
	num_proc = options.num_proc
	bowtie2_build = options.bowtie2_build
	ref_db_fh = tarfile.open(db_tarfile, 'r')
	# make a list of ref ffn files to extract
	files_to_extract = []
	unclassified_sp = {}
	species_ref = {}
	for sp in species_list:
		if sp.count('unclassified') == 1:
			genus = re.search('(.+)_unclassified$', sp).group(1)
			if genus not in unclassified_sp:
				unclassified_sp[genus] = []
	for sp in species_list:
		if sp.count('unclassified') == 1: continue
		genus = sp.split('_')[0]
		if genus in unclassified_sp: unclassified_sp.pop(genus, None)
		species_ref[sp] = None
			
	for info in ref_db_fh:
		if info.name.count('ffn.bz2') == 0:
			continue
		gn, sp = re.search('g__(.+)\.s__(.+)\.centroids\.ffn', info.name).group(1, 2)
		if gn in unclassified_sp:
			unclassified_sp[gn].append(info.name)
		elif sp in species_ref:
			species_ref[sp] = info.name
	
	for gn in unclassified_sp:
		for file in unclassified_sp[gn]:
			files_to_extract.append(file)
	for sp in species_ref:
		if species_ref[sp] == None: continue
		files_to_extract.append(species_ref[sp])
	
	temp_dir = options.outdir + '/sams/'
	if not os.path.exists(temp_dir): os.mkdir(temp_dir)
	bz2_dir = temp_dir + '/db_bz2/'
	if not os.path.exists(bz2_dir): os.mkdir(bz2_dir)
	
	# extract the files using tarball index file
	extract_tarball(db_tarfile, db_tarindex, bz2_dir, files_to_extract)
	
	ref_file = temp_dir + '/merged_ref.ffn'
	ref_prefix = temp_dir + '/merged_ref'
	merged_ref = open(ref_file, 'w')
	
	for ffn in glob.glob(bz2_dir +'/*.ffn.bz2'):
		sp = re.search('(s__.+)\.centroids', ffn).group(1)
		for line in bz2.BZ2File(ffn, 'r'):
			if line[0] == '>':
				line = line[:-1] + '|' + sp +'\n'
			merged_ref.write(line)
	merged_ref.close()
	
	# bowtie2_build index merged ffn, run bowtie2 on input files against index
	bowtie2_build_cmd = [bowtie2_build, ref_file, ref_prefix]
	db_process = Popen(bowtie2_build_cmd, stdout = PIPE, stderr = PIPE)
	p_out, p_err = db_process.communicate()
	
	shutil.rmtree(bz2_dir)

def align_sp(args):
	species_list, db_tarfile, db_tarindex, infiles, sample_name, options = args
	infile1, infile2 = infiles
	bowtie2 = options.bowtie2
	samtools = options.samtools
	num_proc = options.num_proc
	ref_file = options.outdir + '/merged_ref.ffn'
	db_prefix = options.outdir + '/sams/merged_ref'
	sam_dir = options.outdir + '/sams/'
	if not os.path.exists(sam_dir):
		os.mkdir(sam_dir)
	samfiles = []
	
	bowtie2_logfh = open(sam_dir + '/bowtie2.' + str(sample_name) + '.log', 'a')
	bam_prefix = sam_dir + '/reads_vs_ref.'+ str(sample_name)
	bamfile = bam_prefix + '.bam'
	bowtie2_cmd = [bowtie2, '-p', str(num_proc), '-x', db_prefix, '-1', infile1, '-2', infile2]
	if not infile2:
		bowtie2_cmd = [bowtie2, '-p', str(num_proc), '-x', db_prefix, infile1]
	print bowtie2_cmd
	samtools_view = [samtools, 'view', '-bhS', '-']
	samtools_sort = [samtools, 'sort', '-', bam_prefix]
	samtools_index = [samtools, 'index', bamfile]
	p1 = Popen(bowtie2_cmd, stdout = PIPE, stderr = bowtie2_logfh)
	p2 = Popen(samtools_view, stdin = p1.stdout, stdout = PIPE, stderr = bowtie2_logfh)
	p3 = Popen(samtools_sort, stdin = p2.stdout, stdout = PIPE, stderr = bowtie2_logfh)
	p1.stdout.close()
	p2.stdout.close()
	output, err = p3.communicate()
	samtools_index = [samtools, 'index', bamfile]
	p1 = Popen(samtools_index, stderr = bowtie2_logfh, stdout = bowtie2_logfh)
	output, err = p1.communicate()
	bowtie2_logfh.close()
	# record finished samples
	samlog = sam_dir + 'sam.log'
	if os.path.exists(samlog): logfh = open(samlog, 'a')
	else: logfh = open(samlog, 'w')
	logfh.write('%s\n' % sample_name)
	logfh.close()
	
# End of align_sp

def interpret_pstring(ref_base, x):
	res = {'start':0, 'end':0, 'insert':0, 'del':0, 'bases': [0,0,0,0]}
	if x == '': return res
	
	res['end'] = x.count('$')
	x = re.sub(r'\$', '', x)
	res['start'] = x.count('^')
	x = re.sub(r'\^.{1}','',x)
	
	for a in re.findall(r'\+\d+', x):
		res['insert'] += 1
		c = int(a.replace('+', ''))
		x=re.sub(r'\+\d+\w{%i}' % c, '', x)
	for a in re.findall(r'\-\d+', x):
		res['del'] += 1
		c = int(a.replace('-', ''))
		x = re.sub(r'\-\d+\w{%i}' % c, '', x)
		
	res['bases'][base_index(ref_base)] += x.count(',') + x.count('.')
	res['bases'][0] += x.count('A') + x.count('a')
	res['bases'][1] += x.count('C') + x.count('c')
	res['bases'][2] += x.count('G') + x.count('g')
	res['bases'][3] += x.count('T') + x.count('t')
	
	return res

def locateSNPs(args):
	species, pileup_file, ref_seq, snp_file, options = args
	from copy import copy
	# dict for gene length
	gene_lengths = {}
	for pid in ref_seq: gene_lengths[pid] = len(ref_seq[pid][1])
	
	if os.stat(pileup_file).st_size == 0:
		return 0
	min_cov = options.min_cov
	minor_cov = 2
	
	# 0. read the pileup file, record candidates from ref_seq
	candidates = {}
	gene_covs = {}
	positions = {}
	ifh = open(pileup_file, 'r')
	while 1:
		line = ifh.readline().rstrip('\n')
		if not line: break
		pid, pos, ref_base, cstring, pstring = line.split('\t')
		if pid not in gene_covs: gene_covs[pid] = []
		pos = int(pos)
		pileups = [interpret_pstring(ref_base, p) for p in pstring.split('|')]
		covs = [int(c) for c in cstring.split(',')]
		gene_covs[pid].append(covs)
		if pid not in positions: positions[pid] = []
		positions[pid].append(pos)
		count = [0,0,0,0]
		for p in pileups:
			for bindex in range(4): count[bindex] += p['bases'][bindex]
		sorted_base_count = sorted(zip(['A','C','G','T'], count), key = lambda x: x[1], reverse = True)
		if sorted_base_count[1][1] > minor_cov:
			if pid not in candidates: candidates[pid] = []
			candidates[pid].append((pos, pileups))
	ifh.close()
	
	# mask pids and filter SNPs
	# transform covs into np.array
	#for pid in gene_covs: gene_covs[pid] = np.array(gene_covs[pid])
	# 1. mask genes with covered perc < 0.75
	#    for those who pass, estimate per gene median coverage across length, define upper/lower bound
	masked_genes = {}
	#gene_median_cov = []
	#for pid in gene_covs:
	#	if float(gene_covs[pid].shape[0])/gene_lengths[pid] < 0.75:
	#		masked_genes[pid] = 1 # covered perc cal.
	#		continue
	#	gene_median_cov.append(np.percentile(gene_covs[pid], 50, axis = 0))
		#print pid, gene_covs[pid].shape[0], gene_lengths[pid], np.percentile(gene_covs[pid], 50, axis = 0)
	#try:
	#	cov_lower = np.percentile(np.array(gene_median_cov).sum(axis=1), 25)
	#	cov_upper = np.percentile(np.array(gene_median_cov).sum(axis=1), 75)
		#print cov_lower, cov_upper		
	#except:
	#	return 0
	# 2. per gene, mask positions outside of the boundary by cov_lower/upper
	#    if masked_positions > 50% of gene length, mask the gene
	masked_positions = {}
	#for pid in gene_covs:
	#	median = np.percentile(gene_covs[pid].sum(axis=1), 50, axis = 0)
		#print pid, median
	#	if median < cov_lower or median > cov_upper: 
	#		masked_genes[pid] = 1
	#		continue
	#	site_sum = gene_covs[pid].sum(axis=1)
	#	gene_lower = np.percentile(site_sum, 25)
	#	gene_upper = np.percentile(site_sum, 75)
	#	gene_masked_positions = []
	#	num_masked = 0
	#	for index in range(site_sum.shape[0]):
	#		if site_sum[index] > gene_upper or site_sum[index] < gene_lower:
	#			if pid not in masked_positions: masked_positions[pid] = {}
	#			masked_positions[pid][positions[pid][index]] = 1
		#print len(masked_positions[pid]), len(positions[pid])
		#masked_rate = float(M)/L
	#	masked_rate = float(len(masked_positions[pid])) / len(positions[pid])
		#print median, site_sum, gene_lower, gene_upper, num_masked, masked_rate
	#	if masked_rate > 0.5: masked_genes[pid] = 1
	#print masked_genes
	# 3. filter genes with >10% SNPs
	#for pid in candidates:
	#	snp_rate = float(len(candidates[pid]))/gene_lengths[pid]
	#	if snp_rate > 0.1: masked_genes[pid] = 1
	
	# 4. instantiate the SNPs	
	SNPs = []
	for pid in candidates:
		if pid in masked_genes: continue
		for pos, pileups in candidates[pid]:
			if pid in masked_positions and pos in masked_positions[pid]: continue
			S = SNP()
			S.gene = pid
			S.loc = pos
			counts = []
			for p in pileups: counts.append(p['bases'])
			S.counts = counts
			S.counts2freqs(minor_cov)
			S.sort_freqs()
			flag = 0  # flag out SNPs with >2 alleles
			for f in S.sorted_freqs:
				if f[2] > 0:
					flag = 1
					break
			if flag == 0: SNPs.append(S)
			#print flag, pid, pos, pileups, S.sorted_freqs
	
	cPickle.dump(SNPs, open(snp_file, 'wb'))
	if not options.quiet: sys.stdout.write('    [%s] Finished SNPs processing\n' % species)
	return 0


def clean_tempdir(outdir):
	shutil.rmtree(temp_dir + '/db_bz2')
# End of clean_tempdir

############################# SNP-FLOW ########################
def SNP_flow(G, node_layers):
	sources, targets = node_layers[0], node_layers[-1]
	models = []
	# trivial case if t_depth == 1
	if len(node_layers) == 1:
		M = Model()
		for nodeID in sources:
			M.profile.append(G.node[nodeID]['perc'])
			M.strains.append(G.node[nodeID]['seq'])
			if len(M.locs) == 0: M.locs = G.node[nodeID]['locs']
		models.append(M)
		return models
	
	# non-trivial case when t_depth > 1
	# 1. get the locs
	locs = []
	for nodes in node_layers:
		nodeID = nodes[0]
		locs += G.node[nodeID]['locs']
	# 2. get all the combo of source/target 
	combos = list(product(*[sources, targets]))
	all_paths = {}
	for source, target in combos:
		simple_paths = list(nx.all_simple_paths(G, source, target))
		for path in simple_paths: all_paths[tuple(path)] = 1
	
	# 3. generate all possible combinations of path/rel_ab
	# init path and capacity
	profiles = {}
	paths = all_paths.keys()
	# generate profiles
	for index in range(len(paths)):
		node_capacity = {}
		for node, node_data in G.nodes(data=True): node_capacity[node] = node_data['perc']
		profile = []
		while available_path_indices(paths, node_capacity) != []:
			curr_path = paths[index]
			curr_ab = path_flow(curr_path, node_capacity)
			if curr_ab <= 0: break
			profile.append((curr_ab, curr_path))
			update_capacity(paths, node_capacity, index)
			index = update_index(paths, node_capacity, index)
		ukey = unique_key(profile)
		if ukey not in profiles:
			profiles[ukey] = profile
	
	models = []
	for x in profiles.values():
		if x == []: continue
		M = Model()
		M.locs = locs
		# get profile and use profile to infer ML strains
		M.profile = map(itemgetter(0), x)
		#M.strains = M.infer_strains(profile, selected_SNPs)
		
		paths = map(itemgetter(1), x)
		for path in paths:
			seq = ''
			for nodeID in path: seq += G.node[nodeID]['seq']
			M.strains.append(seq)
		models.append(M)
	return models

def available_path_indices(paths, node_capacity):
	indices = []
	for i, p in enumerate(paths):
		if path_flow(p, node_capacity) > 0: indices.append(i)
	return indices

def update_capacity(paths, node_capacity, index):
	flow = path_flow(paths[index], node_capacity)
	for node in paths[index]: node_capacity[node] -= flow
	
def path_flow(path, node_capacity):
	flows = []
	for node in path:
		flows.append(node_capacity[node])
	if sorted(flows)[0] <= 0: return -1
	return sorted(flows)[0]
	
def update_index(paths, node_capacity, index):
	if index < len(paths) - 1:
		for i, p in enumerate(paths[index+1:]):
			if path_flow(p, node_capacity) > 0: return index+i+1
	for i, p in enumerate(paths[:index]):
		if path_flow(p, node_capacity) > 0: return i
	return index

def unique_key(combo):
	subkey = []
	for x in combo:
		s = []
		for y in sorted(x[-1]): s.append(str(y))
		subkey.append(':'.join(s))
	return '&'.join(sorted(subkey))


############################## SNP CLUSTERING ##############################
def SNPdist(a, b):
	N = len(a)/4
	X = np.array(a).reshape((N, 4))
	Y = np.array(b).reshape((N, 4))
	distances = []
	for i in range(N): # pivot around a certain sample
		# sort columns by the i-th row (sample)
		sampleX = X[:, X[i, :].argsort()]
		sampleY = Y[:, Y[i, :].argsort()]
		d = (sampleX - sampleY)**2
		d = np.sqrt(d.sum(axis=1)).sum()/N
		distances.append(d)
	return sorted(distances)[0]
	
def pick_SNP_medoid(D, indices):  # D is the flat form of all-vs-all distance
	N=condense_matrix_dimen(D)
	#print len(D), N
	distances = {}
	for i in indices:
		if i not in distances: distances[i] = 0
		for j in indices:
			if i == j: continue
			distances[i] += condense_matrix_dist(D, N, i, j)
	#print distances
	medoid_index = sorted(distances.iteritems(), key = lambda x: x[1])[0][0]
	return medoid_index
	
def rep_SNP_cluster(SNPs, sample_index):
	prime_seq = ''
	minor_seq = ''
	locs = []
	X = []
	for S in SNPs:
		freqs = S.sorted_freqs[sample_index]
		if min(freqs) < 0: freqs = [0., 0., 0., 0.] # reset to 0. if -1
		X.append(freqs)
		prime_seq += S.sorted_bases[sample_index][0]
		minor_seq += S.sorted_bases[sample_index][1]
		locs.append([S.gene, S.loc])
	F = np.array(X).mean(axis = 0)
	prime_perc = F[0]
	minor_perc = F[1]
	return locs, (prime_perc, minor_perc), (prime_seq, minor_seq)

def clusterSNPs(SNPs, masked_samples, cluster_perc_thr = 0.05, cluster_size_thr = 10):
	# hierarchical clustering first
	# sort every position allele frequency in descending order first
	# decide max_layers with depth
	
	C = np.array([a.counts for a in SNPs])
	avg_C = C.sum(axis=2).mean()
	max_layers = min(10, int(avg_C)/5)
	#print C
	#print avg_C
	num_samples = len(SNPs[0].counts)
	L = len(SNPs)  # L is the number of total SNPs
	X = [] # pan frequencies
	
	for S in SNPs:
		F = []
		for i in range(num_samples):
			if i in masked_samples: continue  # ignore the samples without sufficient coverage
			F += S.sorted_freqs[i]
		#print F
		X.append(F)
	X = np.array(X)
	#print X.shape
	hclust = hierarchy.linkage(X, method = 'complete', metric = 'euclidean')
	#print hclust
	np.clip(hclust, 0., 1e10, hclust)
	# cut the tree at diff depth, ranging from 0.0 to 1.0 with 0.01 step
	prev_fclust = []
	fclust = []
	fclusts = []
	for cluster_tol in np.arange(0., 2, 0.01):
		fclust = hierarchy.fcluster(hclust, t = cluster_tol, criterion='inconsistent')
		if not np.array_equal(fclust, prev_fclust):
			prev_fclust = fclust
			fclusts.append((cluster_tol, fclust))
			#print fclust
		if max(fclust) == 1: break
	
	D = pdist(X)		
	SNP_clusters = []  # dimension: num_depths * num_clusters
	for tol, fclust in fclusts:
		num_clusters = max(fclust)
		# put different SNPs into different clusters
		clusters = []
		for i in range(num_clusters): clusters.append([])
		for SNP_index, clst_ID in enumerate(fclust): # i is the cluster index
			clusters[clst_ID-1].append(SNP_index)
		#sorted_clustrs = sorted(clusters, key = lambda a: len(a), reverse = True)
		# filter clusters using their size
		# rank clusters using size
		filtered_clusters = []
		cluster_medoids = []
		for cluster in clusters:
			#print len(cluster)
			if float(len(cluster))/L < cluster_perc_thr or len(cluster) < cluster_size_thr: continue
			filtered_clusters.append([])
			for SNP_index in cluster: filtered_clusters[-1].append(SNPs[SNP_index])
			# pick the medoid for the cluster
			medoid_index = pick_SNP_medoid(D, cluster)
			cluster_medoids.append(medoid_index)
		#print tol, len(cluster_medoids)
		if len(cluster_medoids) > 0 and len(cluster_medoids) <= max_layers:
			SNP_clusters.append((tol, cluster_medoids, filtered_clusters))
	
	return SNP_clusters

################################# STRAIN CLUSTERING ############################
def strain_dist(x, y):
	concordant = 0
	disconcordant = 0
	for a, b in zip(list(x), list(y)):
		if a == '-' or b == '-': continue
		if a == b: concordant += 1
		else: disconcordant += 1
	if (concordant+disconcordant) == 0: return 1.
	dist = float(disconcordant)/(concordant+disconcordant)
	return dist

def select_rep_geno(genos, clusters):
	D = []
	for x in clusters:
		D.append((genos[x], genos[x].count('-')))
	return sorted(D, key = lambda a: a[1])[0][0]

def condense_matrix_dimen(X):
	try:
		N = int((np.sqrt(1+8*len(X))+1)/2)
		return N
	except:
		return -1

def condense_matrix_dist(X, N, i, j):
	if i == j: return 0
	if i > j: i, j = j, i
	condense_index = (2*N-i-2)*(i+1)/2 - N + j
	return X[condense_index]

def membership(X, N, strainID, medoids):
	distances = []
	for clusterIndex, medoid in enumerate(medoids):
		dist = condense_matrix_dist(X, N, strainID, medoid)
		distances.append((clusterIndex+1, dist))
	strain_membership = sorted(distances, key = lambda x: x[1])[0][0]
	return strain_membership

def memberships_comp(x, y):
	for a, b in zip(x, y):
		if a != b:
			return True
	return False

def pick_medoids(X, memberships):
	clusters = []
	N = condense_matrix_dimen(X)
	for i in range(max(memberships)): clusters.append([])
	for strainID, clusterID in enumerate(memberships):
		clusters[clusterID-1].append(strainID)
	medoids = []
	for cluster in clusters:
		if len(cluster) == 0: continue
		dist = {}
		for x in cluster:
			dist[x] = 0
			for y in cluster:
				dist[x] += condense_matrix_dist(X, N, x, y)
		medoid = sorted(dist.iteritems(), key = lambda a: a[1])[0][0]
		medoids.append(medoid)
	return medoids

def sq_error(X, clusters):
	N = condense_matrix_dimen(X)
	intra_dist = []
	inter_dist = []
	for cluster in clusters:
		for i in cluster:
			for j in cluster:
				intra_dist.append(condense_matrix_dist(X, N, i, j))
	for clusterIndex1, cluster1 in enumerate(clusters):
		for clusterIndex2, cluster2 in enumerate(clusters):
			if clusterIndex1 >= clusterIndex2: continue
			for i in cluster1:
				for j in cluster2:
					inter_dist.append(condense_matrix_dist(X, N, i, j))
	try:
		avg_intra = sum(intra_dist)/len(intra_dist)
		avg_inter = sum(inter_dist)/len(inter_dist)
	except:
		return -1
	return avg_intra/avg_inter


def prefilter_strains(profiles, locs, SNPs_dict, all_strains, strain_clusters):
	strains = []
	numSample, numStr = profiles.shape
	for i in range(numStr): strains.append(['-']*len(locs))
	ps = profiles.tolist()
	
	# fill the positions with certain bases
	for strain_index, str_index_cluster in enumerate(strain_clusters):
		for locale_index in range(len(locs)):
			bases = []
			for str_index in str_index_cluster: bases.append(all_strains[str_index][locale_index])
			if len(set(bases)) == 1 and bases.count('-') == 0:
				strains[strain_index][locale_index] = bases[0]
				
	init_locs = []
	strain_strings = []
	for i in range(len(strain_clusters)): strain_strings.append('')
	for i in range(len(locs)):
		locale_bases = [strains[j][i] for j in range(numStr)]
		if locale_bases.count('-') > 0: continue # trivial case
		loc = locs[i]
		init_locs.append(loc)
		for str_index, x in enumerate(locale_bases):
			strain_strings[str_index]+=x
	
	if len(init_locs) < 0.2*len(locs):
		num_loc = int(0.2*len(locs) - len(init_locs))
		k = 0
		for i in range(len(locs)):
			locale_bases = [strains[j][i] for j in range(numStr)]
			if locale_bases.count('-') == 0: continue # trivial case
			# if there are uncertain bases, use inference
			blank_count = locale_bases.count('-')
			min_err = 100000
			prime_bases = []
			loc = locs[i]
			pos_freqs = np.array(SNPs_dict[tuple(loc)].freqs)
			X = np.array(SNPs_dict[tuple(loc)].counts)
			base_cols = X.sum(axis=0)
			base_indexes = np.nonzero(base_cols)[0].tolist()
			prime_bases = []	
			for proposed_bases in list(product(base_indexes, repeat=blank_count)):
				strain_bases = []
				blank_indices = []
				for index, b in enumerate(locale_bases):
					if b == '-':
						strain_bases.append('-')
						blank_indices.append(index)
					else: strain_bases.append(base_index(b))
				for index, base_ind in zip(blank_indices, proposed_bases):
					strain_bases[index] = base_ind
				err = fit_bases(pos_freqs, strain_bases, ps)
				if err < min_err:
					prime_bases = strain_bases
					min_err = err
			for str_index, base_ind in enumerate(prime_bases):
				strains[str_index][i] = ['A','C','G','T'][base_ind]
			locale_bases = [strains[j][i] for j in range(numStr)]
			init_locs.append(loc)
			for str_index, x in enumerate(locale_bases): strain_strings[str_index]+=x
			k+=1
			if k >= num_loc: break
		
	return strain_strings, init_locs

def infer_base_RSS(cohort_locs, SNPs_dict, masked_samples):
	err_sum = 0
	err_num = 1
	for loc in cohort_locs:
		freqs = SNPs_dict[tuple(loc)].freqs
		err_num += len(freqs)
		min_err = 1e20
		for i in range(4):
			F = [0,0,0,0]
			F[i] = 1
			err = 0
			for freq in freqs:
				if sum(freq) <= 0: continue
				err += freq_JSD(freq, F)
			if err < min_err: min_err = err
		err_sum += min_err
	#return err_sum/err_num
	return err_sum/(len(cohort_locs))
	
def infer_strains(profiles, locs, SNPs_dict, all_strains, strain_clusters):
	strains = []
	numSample, numStr = profiles.shape
	for i in range(numStr): strains.append(['-']*len(locs))
	ps = profiles.tolist()
	
	# fill the positions with certain bases
	for strain_index, str_index_cluster in enumerate(strain_clusters):
		for locale_index in range(len(locs)):
			bases = []
			for str_index in str_index_cluster: bases.append(all_strains[str_index][locale_index])
			if len(set(bases)) == 1 and bases.count('-') == 0:
				strains[strain_index][locale_index] = bases[0]
	
	for i in range(len(locs)):
		locale_bases = [strains[j][i] for j in range(numStr)]
		if locale_bases.count('-') == 0: continue # trivial case
		# if there are uncertain bases, use inference
		blank_count = locale_bases.count('-')
		min_err = 100000
		prime_bases = []
		loc = locs[i]
		pos_freqs = np.array(SNPs_dict[tuple(loc)].freqs)
		X = np.array(SNPs_dict[tuple(loc)].counts)
		base_cols = X.sum(axis=0)
		base_indexes = np.nonzero(base_cols)[0].tolist()
		prime_bases = []
		
		for proposed_bases in list(product(base_indexes, repeat=blank_count)):
			strain_bases = []
			blank_indices = []
			for index, b in enumerate(locale_bases):
				if b == '-':
					strain_bases.append('-')
					blank_indices.append(index)
				else: strain_bases.append(base_index(b))
			for index, base_ind in zip(blank_indices, proposed_bases):
				strain_bases[index] = base_ind
			err = fit_bases(pos_freqs, strain_bases, ps)
			if err < min_err:
				prime_bases = strain_bases
				min_err = err
		for str_index, base_ind in enumerate(prime_bases):
			strains[str_index][i] = ['A','C','G','T'][base_ind]
		locale_bases = [strains[j][i] for j in range(numStr)]
	strain_strings = [''.join(x) for x in strains]
	
	return strain_strings

def strain_Kmeans(X, K):
	max_iter = 1000
	N = condense_matrix_dimen(X)  # N is the number of strains
	# start with K randomly drawn medoids
	iterID = 0
	changed = True
	memberships = []
	for i in range(N): memberships.append(random.randint(1, K))  # randomly assign memberships
	medoids = pick_medoids(X, memberships)  # init medoids
	
	# iterative block to determine memberships of each strain
	while iterID < max_iter and changed == True:
		iterID += 1
		# update memberships
		new_memberships = []
		for strainID in range(N):
			strain_membership = membership(X, N, strainID, medoids)
			new_memberships.append(strain_membership)
		changed = memberships_comp(memberships, new_memberships)
		# update clusters
		memberships = new_memberships
		medoids = pick_medoids(X, memberships)
	
	# generate final clusters
	clusters = []
	for i in range(max(memberships)): clusters.append([])
	for strainID, clusterID in enumerate(memberships):
		clusters[clusterID-1].append(strainID)
	err = sq_error(X, clusters)
	
	return clusters, err

def evaluate_cohort_model(strains, profiles, locs, SNPs_dict, masked_samples):
	err = 0
	j = 0
	for sample_index, profile in enumerate(profiles.tolist()):
		if sample_index in masked_samples or sum(profile) == 0: continue
		j += 1
		counts = []
		for i in range(len(locs)):
			S = SNPs_dict[tuple(locs[i])]
			counts.append(S.counts[sample_index])
		M, O, F = to_matrices(strains, counts)
		err += fit_err(np.array(profile), M, F)
	RSS = err/j
	return RSS	

################################### ML fitting ##############################

def JS_entropy(x):
	# jsd is 1/2(KL(a||m)+KL(b||m))
	m = [(i+j)/2 for i, j in zip(x[:4], x[4:])]
	jsd = 0.5*(entropy(x[:4], m) + entropy(x[4:], m))
	return jsd	

def freq_JSD(x, y):
	m = [(i+j)/2 for i, j in zip(x, y)]
	jsd = 0.5*(entropy(x, m) + entropy(y, m))
	return jsd	

def simplex_uniform(N):
	X = uniform(0., 1., N)
	return (-np.log(X))/np.sum(-np.log(X))

# given A, and sigma,  
def gaussian_pick_simplex(A, sigma = SIGMA):
	exp_phis = []
	for a in A.tolist():
		exp_phis.append(np.exp(normal(a, sigma)))
	return np.array(exp_phis)/sum(exp_phis)
	
def fit_err(A, M, F):
	N, L = M.shape[1:]
	O_prime = np.zeros((4, L))
	for i in range(4): O_prime[i] = A.dot(M[i,])
	D = np.zeros((8, L))
	for i in range(L):
		D[4:, i] = O_prime[:, i]
		D[:4, i] = F[:, i]
	# err as JSD
	e = np.apply_along_axis(JS_entropy, 0, D)
	err = np.ma.masked_invalid(e).mean()
	return err

def to_matrices(strs, counts):
	N = len(strs)   # number of strains
	L = len(strs[0])  # number of SNPs
	M = np.zeros((4, N, L)) 
	O = np.zeros((4, L))
	
	for j in range(N):
		for k in range(L):
			i = base_index(strs[j][k])
			M[i, j, k] = 1
	
	for i in range(L):
		a, c, g, t = counts[i]
		O[0, i] = a
		O[1, i] = c
		O[2, i] = g
		O[3, i] = t
	x = O.sum(axis = 0) # coverage along the SNPs
	F = np.zeros((4, L))
	for j in range(L):
		if x[j] == 0: continue
		for i in range(4):
			F[i, j] = O[i, j]/x[j]
	return M, O, F
	
def Metropolis_Hastings(strains, counts, max_iter):
	numStr = len(strains)
	M, O, F = to_matrices(strains, counts)
	profiles = []
	for i in range(max_iter):
		p = simplex_uniform(numStr)  # randomly init profile
		err = fit_err(p, M, F)
		profiles.append((p, err))
	profile = sorted(profiles, key = lambda x: x[1])[0][0]
	return profile
	
################################### STRAIN FITTING #############################
def fit_strains(options, samples):
	N = len(samples)
	species_ab, species_cov = cPickle.load(open(options.outdir + '/species.list', 'rb'))
	species_list = {}  # the ones qualified for further analysis
	for sp in species_cov:
		if sp.count('unclassified'): continue
		if max(species_cov[sp]) >= options.min_cov: species_list[sp] = species_cov[sp]
	
	# load ref sequences
	ref_seqs = {}
	merged_ref = options.outdir + '/merged_ref.ffn'
	for record in SeqIO.parse(merged_ref, 'fasta'):
		species = record.name.split('|')[-1].replace('s__','')
		pid = record.name.split('|')[-2]
		if species not in ref_seqs: ref_seqs[species] = {}
		ref_seqs[species][pid] = record.seq
	
	fit_cmds = []
	species_SNPs_num = {}
	snpflow_dir = options.outdir+'/snpflows/'
	for species in species_list:
		sp_cov = species_list[species]
		species_file = options.outdir + '/pileups/' + species+'.pileups'
		if (not os.path.exists(species_file)) or os.stat(species_file).st_size == 0: continue
		model_file = snpflow_dir + species + '.model'
		snp_file = snpflow_dir + species + '.snps'
		if os.path.exists(model_file) or (not os.path.exists(snp_file)): continue
		try:SNPs = cPickle.load(open(snp_file, 'rb'))
		except: sys.stderr.write('  [WARNING]: Error in detecting SNPs for %s. Skipping...\n' % species)		
		if len(SNPs) < 10:
			sys.stderr.write('  [WARNING]: There are only %i informative SNPs detected for %s.\n' % (len(SNPs), species))
			sys.stderr.write('             You need 20 or more SNPs to confidently infer strains. Skipping...\n')
			continue
		cmd = [species, sp_cov, SNPs, ref_seqs[species], options]
		species_SNPs_num[species] = len(SNPs)
		fit_cmds.append(cmd)

	if not options.quiet:
		sys.stdout.write('  ##### Number of informative SNPs for species #####\n')
		for species in species_SNPs_num: sys.stdout.write('    %s\t%i\n' % (species, species_SNPs_num[species]))
		sys.stdout.write('  ##################################################\n\n')
	
	if not options.quiet:
		sys.stdout.write("  Fitting strains in different species ...\n")
	# mp mode
	"""
	if len(fit_cmds) > 0:
		pool = mp.Pool(options.num_proc)
		pool.map_async(fit_snpflow, fit_cmds)
		pool.close()
		pool.join()
	"""
	for cmd in fit_cmds: fit_snpflow(cmd)  # single thread mode for testing
	
	if not options.quiet:
		sys.stdout.write('  Done.\n')
		
# End of fit_strains

def ML_strain(SNPs, profile, sampleIndex):
	strains = []
	error = 0
	N = len(profile)
	for i in range(N): strains.append('')
	for SNP_index, S in enumerate(SNPs):
		freqs = S.freqs[sampleIndex]
		if freqs[0] == -1:
			x = []
			for i in range(N): x.append('-')
			strains.append(x)	
			continue
		sorted_bases = S.sorted_bases[sampleIndex]
		sorted_freqs = S.sorted_freqs[sampleIndex]
		bases = {}
		for b, f in zip(sorted_bases, sorted_freqs): bases[b] = f
		# all combinations of bases
		min_err = 1e5
		optimal_genotype = None
		for genotype in list(product(bases.keys(), repeat = N)):
			X = {}
			for g, f in zip(genotype, profile):
				if g not in X: X[g] = f
				else: X[g] += f
			exp_freqs = []
			for b in ['A','C','G','T']:
				if b not in X: exp_freqs.append(0)
				else: exp_freqs.append(X[b])
			err = euclidean(freqs, exp_freqs)
			if err < min_err:
				min_err = err
				optimal_genotype = genotype
		error += min_err
		for i in range(N):
			strains[i] += optimal_genotype[i]
	return strains

def gen_SNP_flow(SNP_clusters, N, masked_samples):
	# SNP-flow estimates possible strain profiles for each sample
	strain_models = []
	for sample_index in range(N):
		if sample_index in masked_samples: continue
		sample_models = []
		for tol, medoids, clusters in SNP_clusters:
			sample_graph = nx.DiGraph()
			sample_node_layers = []
			nodeID = 0
			# initiate the graphs and node_layers
			# 1. add nodes and attributes
			for medoid, cluster in zip(medoids, clusters):
				locs, freqs, seqs = rep_SNP_cluster(cluster, sample_index)
				if freqs[1] <= 0: break
				sample_node_layers.append([])
				# normalize freqs
				major_freq = freqs[0]/sum(freqs[:2])
				minor_freq = 1 - major_freq
				nodeID += 1
				sample_node_layers[-1].append(nodeID)
				#print sample_index, medoid, len(clusters), len(cluster), major_freq, minor_freq
				sample_graph.add_node(nodeID, {'perc': major_freq, 'seq': seqs[0], 'locs': locs})
				nodeID += 1
				sample_node_layers[-1].append(nodeID)
				sample_graph.add_node(nodeID, {'perc': minor_freq, 'seq': seqs[1], 'locs': locs})
			# 2. add edges to graph
			if len(sample_node_layers) >= 2 and nodeID > 0:
				for node_set1, node_set2 in zip(sample_node_layers[:-1], sample_node_layers[1:]):
					x = [node_set1, node_set2]
					for node1, node2 in list(product(*x)):
						sample_graph.add_edge(node1, node2)
			# 3. append models
			if len(sample_graph.nodes()) > 0:
				models = SNP_flow(sample_graph, sample_node_layers)
				sel_models = []
				for model in models:
					if len(model.profile) == 0: continue
					sel_models.append(model)
				strain_models.append((sample_index, tol, sel_models))
	
	return strain_models


def fit_snpflow(args):
	species, sp_cov, SNPs, ref_seq, options = args
	species_file = options.outdir + '/pileups/' + species+'.pileups'
	sample_model_file = options.outdir + '/snpflows/' + species + '.sample_models'
	model_file = options.outdir + '/snpflows/' + species + '.model'
	np.seterr(all='ignore')
	if not options.quiet:
		sys.stdout.write('  [%s] 1. Constructing strain models using SNP-flow...\n' % species)
	# masking samples where the coverage is not deep enough for inference
	masked_samples = [i for i in range(len(sp_cov)) if sp_cov[i] < options.min_cov]
	
	# sort SNPs with overall counts and select the top 1000, or max num of SNPs
	selected_SNPs = sorted(SNPs, key = lambda x: x.total_count(), reverse = True)[:1000]
	locs = []
	for S in selected_SNPs: locs.append((S.gene, S.loc))
	N = len(selected_SNPs[0].counts)  # N = number of samples
	L = len(selected_SNPs)	# number of SNPs
	#print locs
	if not os.path.exists(sample_model_file):
		# cluster SNPs using all samples
		SNP_clusters = clusterSNPs(selected_SNPs, masked_samples)
		strain_models = gen_SNP_flow(SNP_clusters, N, masked_samples)
		
		if not options.quiet:
			sys.stdout.write('  [%s] 2. %i candidate models constructed.\n' % (species, len(strain_models)))
		if len(strain_models) == 0:
			return 0
	
		# 1. go through each model to evaluate
		strains = {}
		sample_optimal_models = {}
		for sample_index in range(N):
			if sample_index in masked_samples: 
				sample_optimal_models[sample_index]=None
				continue
			selected_models = []
			all_sample_models = {}  # sort all models by number of strains
			for i, tol, models in strain_models:
				if i != sample_index: continue
				for model in models:
					numStrains = len(model.profile)
					if numStrains not in all_sample_models: all_sample_models[numStrains] = []
					all_sample_models[numStrains].append(model)
		
			optimal_model = None
			min_AICc = 1e5
			for numStrains in sorted(all_sample_models.keys()):
				models = all_sample_models[numStrains]
				local_min_AICc = 1e5
				local_optimal_model = None
				for model in random.sample(models, min(32, len(models))):
					with np.errstate(invalid='ignore'): model.evaluate_model(selected_SNPs, sample_index)
					model.calAICc(numStrains, sample_index)
					if model.AICc < local_min_AICc:
						local_min_AICc = model.AICc
						local_optimal_model = model
				if local_min_AICc < min_AICc:
					if len(models) <= 32:
						optimal_model = local_optimal_model
						min_AICc = local_min_AICc
					else:
						for model in models:
							with np.errstate(invalid='ignore'): model.evaluate_model(selected_SNPs, sample_index)
							model.calAICc(numStrains, sample_index)
							if model.AICc < min_AICc:
								optimal_model = model
								min_AICc = model.AICc
			
			try: aligned_strains = optimal_model.align_strains(locs)
			except: optimal_model = None
			sample_optimal_models[sample_index]=optimal_model
		cPickle.dump(sample_optimal_models, open(sample_model_file, 'wb'))
		if not options.quiet:
			sys.stdout.write('  [%s] 3. Done constructing sample specific optimal models.\n' % species)
	
	else:
		if not options.quiet:
			sys.stdout.write('  [%s] 2-3. Loading optimal models for each sample in cohort...\n' % species)
		sample_optimal_models = cPickle.load(open(sample_model_file, 'rb'))
		
	# 2. unify the sample models into cohort model	
	#  computationally expensive if N is large, solve with dynamic programming	
	N = len(selected_SNPs[0].counts)  # N = number of samples
	if N == 1: # trivial case of N == 1:
		model_prime = sample_optimal_models[0]
		cohort_strains = model_prime.strains
		init_profiles = [model_prime.profile]
		cohort_locs = model_prime.aligned_locs
		# init seed_profile
		SNPs_dict = {} # init dict for SNPs keyed by loc
		for S in SNPs: SNPs_dict[(S.gene, S.loc)] = S
		seed_profile = np.zeros((N, len(cohort_strains))) # row as sample and col as strain
		for i in range(N):
			if i in masked_samples: continue
			for j in range(len(cohort_strains)):
				seed_profile[i,j] = init_profiles[i][j]
	else: # non-trivial case
		# 2.1 get all strains
		strain_dict = {}
		cohort_locs = []
		all_strains = []
		for sample_index in range(N):
			if sample_optimal_models[sample_index] == None: continue
			for strain_index, strain in enumerate(sample_optimal_models[sample_index].aligned_strains):
				all_strains.append(strain)
				if cohort_locs == []: cohort_locs = sample_optimal_models[sample_index].aligned_locs
				strain_rel_ab = sample_optimal_models[sample_index].profile[strain_index]
				if strain not in strain_dict: strain_dict[strain] = [(sample_index, strain_index, strain_rel_ab)]
				else: strain_dict[strain].append((sample_index, strain_index, strain_rel_ab))
		
		# auto set upper-bound K
		Kmax = len(all_strains)
		# build the tree using all strains
		X = [] # X is the compressed distance matrix
		for i, strain_a in enumerate(all_strains):
			for j, strain_b in enumerate(all_strains):
				if i >= j: continue
				X.append(strain_dist(strain_a, strain_b))
		Ks_clustering = {}
		for numStr in range(2, Kmax+1):
			# generate combinations
			for r in range(10): # repeat 10 times to mitigate randomness
				strain_clusters, sq_err = strain_Kmeans(X, numStr)
				k = len(strain_clusters)
				if sq_err == -1: continue
				if k not in Ks_clustering: Ks_clustering[k] = [strain_clusters, sq_err]
				elif sq_err < Ks_clustering[k][1]:Ks_clustering[k] = [strain_clusters, sq_err]
		
		# pick the representative model
		optimal_model = []
		optimal_K = 0
		min_AICc = 1000000
		SNPs_dict = {} # init dict for SNPs keyed by loc
		for S in SNPs: SNPs_dict[(S.gene, S.loc)] = S
		
		# get profiles of the cohort, keyed by number of strains
		candidate_profiles_by_numStr = {}
		for numStrains in sorted(Ks_clustering.keys()): # go thru numStr in ascending order
			clusters, err = Ks_clustering[numStrains]
			# get seed profiles in candidate_profiles
			candidate_profiles = np.zeros((N, numStrains))  # each row is a sample
			for rep_strain_index, cluster in enumerate(clusters):
				rep_rel_ab = {}  # records acc. rel_ab for rep strain, keyed by sample index
				for i in cluster:
					strain = all_strains[i]
					for sample_index, strain_index, strain_rel_ab in strain_dict[strain]:
						if sample_index not in rep_rel_ab: rep_rel_ab[sample_index] = 0
						rep_rel_ab[sample_index] += strain_rel_ab
						print strain_rel_ab
				# update candidate_profiles using rep_rel_ab
				for sample_index in rep_rel_ab:
					candidate_profiles[sample_index, rep_strain_index] = rep_rel_ab[sample_index]
			# normalize profiles
			for sampleIndex in range(candidate_profiles.shape[0]):
				if candidate_profiles[sampleIndex,].sum() == 0: continue
				candidate_profiles[sampleIndex,] /= candidate_profiles[sampleIndex,].sum()
			candidate_profiles_by_numStr[numStrains] = [candidate_profiles, clusters]
			print candidate_profiles
			
		# fast pre-filter proposed models
		# get the base RSS
		base_RSS = infer_base_RSS(cohort_locs, SNPs_dict, masked_samples)
		print base_RSS
		RSSs = [(1, base_RSS)]
		for numStrains in sorted(candidate_profiles_by_numStr.keys()):
			candidate_profiles, clusters = candidate_profiles_by_numStr[numStrains]
			# get seed strains using ML inference
			init_strains, init_locs = prefilter_strains(candidate_profiles, cohort_locs, SNPs_dict, all_strains, clusters) # expensive part
			#print init_strains, init_locs
			# evaluate candidate_model using AICc
			n = len(init_locs)
			RSS = evaluate_cohort_model(init_strains, candidate_profiles, init_locs, SNPs_dict, masked_samples)	
			RSSs.append((numStrains, RSS))
		print RSSs
		# get the optimal model
		try:
			optimal_numStrains = 1
			for (x1, y1), (x2, y2) in zip(RSSs[:-1], RSSs[1:]):
				delta = (y1 - y2)/y1
				print delta
				if delta <= 0.1:  # if less than 5% improvement 
					optimal_numStrains = x1
					break
				else: optimal_numStrains = x2
			print optimal_numStrains
			if optimal_numStrains == 1:
				return 0
			seed_profile, optimal_clusters = candidate_profiles_by_numStr[optimal_numStrains]
			cohort_strains = infer_strains(seed_profile, cohort_locs, SNPs_dict, all_strains, optimal_clusters) # expensive part
			#print seed_profile
			#print optimal_clusters
			#print cohort_strains
		except ValueError: return 0
		if not options.quiet: sys.stdout.write('  [%s] 4. Unified sample models into cohort model with %s strains.\n' % (species, len(cohort_strains)))
		
		if not options.quiet: sys.stdout.write('  [%s] 5. All steps successfully completed.\n' % species)
	


# End of fit_snpflow
	
################################### CONSTRUCT WORKFLOW ######################################
def construct(samples, db_tarfile, db_tarindex, options):
	# Run metaphlan profiling and select species that pass options.depth
	sp_list_file = options.outdir + '/species.list'
	
	if os.path.exists(sp_list_file):
		species_ab, species_cov = cPickle.load(open(sp_list_file, 'rb'))
	else:
		species_ab, species_cov = metaphlan2_filter(samples, options)
		cPickle.dump([species_ab, species_cov],  open(options.outdir+'/species.list', 'wb'))
	
	species_list = {}  # the ones qualified for further analysis
	for sp in species_cov:
		if max(species_cov[sp]) >= options.min_cov: species_list[sp] = species_cov[sp]
	
	if not options.quiet:
		sys.stdout.write("There are %i species selected for straining.\n" % len(species_list))
		sys.stdout.write("  Selected species list and sample coverage:\n    # Species\t%s\n" % '\t'.join([x.name for x in samples]))
	for sp in sorted(species_list):
		cov_strings = []
		for i in range(len(samples)):
			cov_strings.append('%.2f' % species_list[sp][i])
		cov_string = '\t'.join(cov_strings)
		sys.stdout.write('    %s: %s\n' % (sp, cov_string))
	if not options.quiet: sys.stdout.write('\n')
	
	# Extract ref sequences from species_list, and run Bowtie2, return samfile and ref db
	# also split BAM files based on the species
	
	# skip aligning if already pileups 
	pileup_dir = options.outdir+'/pileups/'
	if len(glob.glob(pileup_dir+'*.pileups')) == 0 or (not os.path.exists(pileup_dir)):
		samdir = options.outdir + '/sams/'
		samlog = samdir + 'sam.log'	
		finished_bams = []
		if os.path.exists(samdir) and os.path.exists(samlog):
			for line in open(samlog, 'r'):
				finished_bams.append(line[:-1])
			
		if not options.quiet:
			sys.stdout.write("Now Bowtie2 aligning reads to references, and call SNPs\n")
		aln_cmds = []
	
		for sample in samples:
			if sample.name not in finished_bams:
				if sample.fq1 != None: sample_infiles = [sample.fq1, sample.fq2]
				else: sample_infiles = [sample.fq,""]
				aln_cmds.append((species_list, db_tarfile, db_tarindex, sample_infiles, sample.name, options))
			
		if len(aln_cmds) > 0:
			# make refDB first
			make_refDB(species_list.keys(), db_tarfile, db_tarindex, options)
			pool = mp.Pool(options.num_proc)
			pool.map_async(align_sp, aln_cmds)
			pool.close()
			pool.join()		
			
		if not options.quiet:
			sys.stdout.write("Done.\n")
	
	# from bowtie2 alignments, construct SNP graph with ref guidance
	if not options.quiet:
		sys.stdout.write("Now Samtools mpileup to call SNPs, and construct SNP flows.\n")
	snpflow_dir = options.outdir+'/snpflows/'
	construct_snpflow(options, samples, species_list)
	if not options.quiet:
		sys.stdout.write("Done extracting SNPs and SNP flow construction.\n")

	os.system('cp %s %s' % (options.outdir+'/sams/merged_ref.ffn', options.outdir))
	
	return 0

################################### OUTPUT ######################################
def chunks(l, n):
	for i in xrange(0, len(l), n): yield l[i:i+n]

def cluster_bases2infer(candidate_infer_bases, profiles, min_cov):
	rep_counts = []
	cluster_seqs = []
	rep_locs = []
	
	i = 0
	# break into pieces of 2Ks
	base_inference={}
	infer2bases_chunks = list(chunks(candidate_infer_bases, 2000))
	for infer2bases_chunk in infer2bases_chunks:
		data = []
		x_bases=[]
		locs = []
		ref_bases = []
		for pid, pos, ref_base, pos_counts in infer2bases_chunk:
			pos_counts = np.array(pos_counts)
			agg_count = pos_counts.sum(axis = 0)
			sorted_counts = sorted(zip(range(4), agg_count), key = lambda x: x[1], reverse = 1)
			major_base, minor_base = sorted_counts[0][0], sorted_counts[1][0]
			x_bases.append([major_base, minor_base])
			locs.append([pid, pos])
			ref_bases.append(ref_base)
			c = pos_counts.sum(axis = 1)
			c[c==0] = 1 # avoid runtimewarning
			pos_freq=pos_counts.astype(np.float16)/c[:, np.newaxis]
			d = pos_freq[:,major_base].tolist() + pos_freq[:,minor_base].tolist()
			data.append(d)
			
		# cluster data
		data = np.array(data)
		try:
			with np.seterr(all='ignore'): kcentroids, data_labels = vq.kmeans2(data, 50)
		except:
			kcentroids, data_labels = vq.kmeans2(data, 50)
		fclusters = []
		for i in range(max(data_labels)): fclusters.append([])
		for data_index, data_label in enumerate(data_labels):
			fclusters[data_label-1].append(data_index)
			
		pos_counts_rep = []
		for cluster in fclusters:
			if len(cluster) == 0: continue
			rep_index = cluster[0]
			rep_major, rep_minor = x_bases[rep_index]
			rep_pos_count = infer2bases_chunk[cluster[0]][3]
			inferred_bases = infer_base(profiles, rep_pos_count, min_cov)
			if '-' in inferred_bases: # trivial case
				for data_index in cluster:
					pid, pos = locs[data_index]
					ref_base = ref_bases[data_index]
					if pid not in base_inference: base_inference[pid] = []
					base_inference[pid].append([pos, ref_base, inferred_bases])
			else:
				inferred_base_index = [base_index(b) for b in inferred_bases]
				base_marks = []
				for x in inferred_base_index:
					if x == rep_major: base_marks.append(0)
					else: base_marks.append(1)
				for data_index in cluster:
					pid, pos = locs[data_index]
					ref_base = ref_bases[data_index]
					if pid not in base_inference: base_inference[pid] = []
					d_major, d_minor = x_bases[data_index]
					inferred_bases = []
					for m in base_marks:
						if m == 0: inferred_bases.append(index_base(d_major))
						else: inferred_bases.append(index_base(d_minor))
					base_inference[pid].append([pos, ref_base, inferred_bases])
	
	return base_inference

def fit_bases(freqs, bases, profiles):
	err = 0
	for i, profile in enumerate(profiles):
		if sum(profile) == 0: continue
		freq = freqs[i,:]
		if len(freq) != 4 or sum(freq) <= 0: continue
		F = [0,0,0,0]
		for base_index, f in zip(bases, profile):
			F[base_index] += f
		err += freq_JSD(freq, F)
	return err
	
def infer_base(profiles, pos_counts, min_cov):
	numStr = len(profiles[0])
	bases = []
	X = np.array(pos_counts) # rows are samples and cols are bases
	if X.sum(axis=1).max() < min_cov/2:
		for i in range(numStr): bases.append('-')
		return bases
	pos_freqs = np.zeros(X.shape)
	for i in range(X.shape[0]):
		for j in range(X.shape[1]):
			try: pos_freqs[i,j] = float(X[i,j])/X.sum(axis=1)[i]
			except: continue
	# infer pos_freqs
	base_cols = X.sum(axis=0)
	base_indexes = np.nonzero(base_cols)[0].tolist()
	if len(base_indexes) > 3:
		for i in range(numStr): bases.append('-')
		return bases
	min_err = 100000
	prime_bases = []
	for strain_bases in list(product(base_indexes, repeat=numStr)):
		with np.errstate(all='ignore'):
			err = fit_bases(pos_freqs, strain_bases, profiles)
		if err < min_err:
			prime_bases = strain_bases
			min_err = err
	for base_index in prime_bases: bases.append(['A','C','G','T'][base_index])
	return bases


	# cluster candidate_infer_bases first
	base_inference = cluster_bases2infer(candidate_infer_bases, profiles, min_cov)
	# merge inferenced bases to uniGcode
	for pid in base_inference:
		if pid not in uniGcode: uniGcode[pid] = base_inference[pid]
		else: uniGcode[pid]+=base_inference[pid]
		
	# write to uniGcode output file
	ufh = open(uniGcode_file, 'w')
	ufh.write('# *: not covered base; -: uncertain base\n')
	ufh.write('#pid\tposition\tref\t%s\n' % '\t'.join(['str-%i' % (i+1) for i in range(numStr)]))
	for pid in sorted(uniGcode.keys()):
		for pos, ref_base, bases in sorted(uniGcode[pid], key=lambda x: x[0]):
			ufh.write('%s\t%i\t%s\t%s\n' % (pid, pos, ref_base, '\t'.join(bases)))
	ufh.close()
	
	if not quiet: sys.stdout.write('  [%s] 3. All finished.\n' % species)
	return 0
	
	
def output_results(samples, options):
	results_dir = options.outdir + '/results/'
	metaphlan_dir = options.outdir + '/metaphlan/'
	projInfo = cPickle.load(open(options.outdir+'/proj.log', 'rb'))
	
	# read rel_ab from metaphlan files
	rel_ab, sp_cov = cPickle.load(open(options.outdir+'/species.list', 'rb'))
	if not os.path.exists(results_dir): os.mkdir(results_dir)

	if not options.quiet:
		sys.stdout.write('  Writing intra-species and strain relative abundance tables...\n')
	# prepare rel. ab. and intra-species ab
	sp_intra_profile = {}
	sp_overall_profile = {}
	sp_masked_samples = {}
	sp_strains = {}
	
	insuf_sp = {}
	single_str_sp = {}
	unresolved_sp = {}
	for species in sp_cov:
		solution_file = options.outdir + '/snpflows/' + species + '.model'
		snp_file = options.outdir + '/snpflows/' + species + '.snps'
		pileup_file = options.outdir + '/pileups/' + species + '.pileups'
		
		if not os.path.exists(pileup_file) and species.count('unclassified') == 0:
			insuf_sp[species] = rel_ab[species]
			continue
		if not os.path.exists(snp_file) and species.count('unclassified') == 0:
			sp_intra_profile[species] = [1.0]*len(samples)
			sp_overall_profile[species] = rel_ab[species]
			unresolved_sp[species] = 1
			continue
		if species.count('unclassified') == 1:
			unresolved_sp[species] = 1
			continue
		SNPs = cPickle.load(open(snp_file, 'rb'))
		if len(SNPs) < 10:
			sp_intra_profile[species] = [1.0]*len(samples)
			sp_overall_profile[species] = rel_ab[species]
			unresolved_sp[species] = 1
			continue
		if not os.path.exists(solution_file):
			unresolved_sp[species] = 1
			continue
			
		strains, profiles, locs = cPickle.load(open(solution_file, 'rb'))
		sp_strains[species] = [locs, strains, profiles]
		M = len(strains)  # M is number of strains
		species_overall_ab = rel_ab[species]
		
		# mark samples w/o sufficient coverage
		sp_masked_samples[species] = [i for i in range(len(samples)) if sp_cov[species][i] < options.min_cov]
		
		# read relative abundance profiles
		sp_intra_profile[species] = []
		sp_overall_profile[species] = []
		numStr = len(strains)
		for str_index in range(numStr):
			sp_intra_profile[species].append([0.]*len(samples))
			sp_overall_profile[species].append([0.]*len(samples))
		for sample_index, sample_profile in enumerate(profiles):
			for str_index in range(numStr):
				overall_rel_ab = sample_profile[str_index]*species_overall_ab[sample_index]
				sp_overall_profile[species][str_index][sample_index] = overall_rel_ab
				sp_intra_profile[species][str_index][sample_index] = 100*sample_profile[str_index]
	
	# output aggregated intra-species profile and overall strain profile
	aggregated_intra_file = results_dir + 'Intra_sp_rel_ab.profiles'
	aggregated_overall_file = results_dir + 'Overall_rel_ab.profiles'
	intrafh = open(aggregated_intra_file, 'w')
	intrafh.write('# Species\tstrain_ID\tmasked_samples')
	relfh = open(aggregated_overall_file, 'w')
	relfh.write('# Species\tstrain_ID\tmasked_samples')
	for i in range(len(samples)):
		intrafh.write('\t%s' % samples[i].name)
		relfh.write('\t%s' % samples[i].name)
	intrafh.write('\n')
	relfh.write('\n')
	
	# write details about species with resolved strains
	for sp in sp_overall_profile:
		try:
			overall_profile = sp_overall_profile[sp]
			intra_profile = sp_intra_profile[sp]
			masked_samples = sp_masked_samples[sp]
		except KeyError: continue
		numStr = len(overall_profile)
		if len(masked_samples) == 0: mask_str = 'NA'
		else: mask_str = ','.join(['%i' % (x+1) for x in masked_samples])
		for str_index in range(numStr):
			overall_profile_str = '\t'.join(['%.8f' % x for x in overall_profile[str_index]])
			intra_profile_str = '\t'.join(['%.8f' % x for x in intra_profile[str_index]])
			relfh.write('%s\tstr-%i\t%s\t%s\n' % (sp, str_index+1, mask_str, overall_profile_str))
			intrafh.write('%s\tstr-%i\t%s\t%s\n' % (sp, str_index+1, mask_str, intra_profile_str))
	# write details about species with insufficient coverage
	for sp in insuf_sp:
		overall_profile_str = '\t'.join(['%.8f' % x for x in rel_ab[sp]])
		intra_profile_str = '\t'.join(['100.00']*len(samples))
		mask_str = ','.join([str(i+1) for i in range(len(samples))])
		relfh.write('%s\tNA,insufficient\t%s\t%s\n' % (sp, mask_str, overall_profile_str))
		intrafh.write('%s\tNA,insufficient\t%s\t%s\n' % (sp, mask_str, intra_profile_str))
	# write details about species with unresolved strains
	for sp in unresolved_sp:
		overall_profile_str = '\t'.join(['%.8f' % x for x in rel_ab[sp]])
		intra_profile_str = '\t'.join(['100.00']*len(samples))
		mask_str = ','.join([str(i+1) for i in range(len(samples))])
		relfh.write('%s\tNA,unresolved\t%s\t%s\n' % (sp, mask_str, overall_profile_str))
		intrafh.write('%s\tNA,unresolved\t%s\t%s\n' % (sp, mask_str, intra_profile_str))
	intrafh.close()
	relfh.close()
	
	if not options.quiet:
		sys.stdout.write('  Done. \n')
		
	
# End of output_results

################################### MAIN ######################################
def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)

	# Compulsory arguments
	compOptions = OptionGroup(parser, "Compulsory parameters",
						"There options are compulsory, and may be supplied in any order.")
	
	compOptions.add_option("-o", "--outdir", type = "string", metavar = "DIR",
							help = "The output directory of Strain-GeMS.")
	
	compOptions.add_option("-c", "--config", type = "string", metavar = "FILE",
							help = "The configuration file of the Strain-GeMS project.")
	
	parser.add_option_group(compOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")

	optOptions.add_option("-t", "--num_proc", type = "int", default = 1, metavar = 'INT',
							help = "Number of processor for ConStrain to use [default: 1].")
							
	optOptions.add_option("-d", "--ref_db", type = "string", metavar = 'STRING',
							help = "The prefix of species reference. [default: Strain-GeMS/db/ref_db].")
							
	optOptions.add_option("-g", "--gsize_db", type = "string", metavar = 'STRING',
							help = "The directory of species average genome size DB. [default: Strain-GeMS/db/gsize.db].")
												
	optOptions.add_option("--bowtie2", type = "string", default = 'bowtie2', metavar = "STRING",
							help = "Path to bowtie2 binary, specify if not in env path [default: bowtie2].\
							Bowtie2 citation: Langmead B. and Salzberg S., Nat. Methods, 2012.\
							Bowtie2 page: http://bowtie-bio.sourceforge.net/bowtie2")
	
	optOptions.add_option("--bowtie2_build", type = "string", default = 'bowtie2-build', metavar = "STRING",
							help = "Path to bowtie2-build binary, specify if not in env path [default: bowtie2-build].\
							Bowtie2 citation: Langmead B. and Salzberg S., Nat. Methods, 2012.\
							Bowtie2 page: http://bowtie-bio.sourceforge.net/bowtie2")
	
	optOptions.add_option("--samtools", type = "string", default = 'samtools', metavar = "STRING",
							help = "Path to samtools binary, specify if not in env path [default: samtools].\
							Samtools citation: Li H., et al, Bioinformatics, 2009.\
							Samtools webpage: http://samtools.sourceforge.net")
	
	optOptions.add_option("-m", "--metaphlan2", type = "string", default = "metaphlan2.py", metavar = "STRING",
							help = "Path to metaphlan2 script, \"metaphlan2.py\", specify if not in env path [default: metaphlan2.py].\
							MetaPhlAn2 citation: Segata N. et al, Nat. Methods, 2012.\
							MetaPhlAn2 page: https://bitbucket.org/biobakery/metaphlan2")
	
	parser.add_option_group(optOptions)
	
	# Parameters that could fine tune the process
	strOptions = OptionGroup(parser, "Straining parameters",
						"There options are optional, and may be supplied in any order.")
	
	strOptions.add_option("--min_cov", type = "float", default = 10, metavar = 'FLOAT',
							help = "Minimum coverage of a species in a sample to be considered [default: 10, range: 5+].")
	
	parser.add_option_group(strOptions)

	# runtime settings that could affect the file saving and message printing
	runtimeSettings = OptionGroup(parser, "Runtime settings",
						"There options are optional, and may be supplied in any order.")
						
	runtimeSettings.add_option("-q", "--quiet", default = False, action = "store_true",
								help = "Suppress printing detailed runtime information, \
								only important warnings/errors will show [default: False].")

	parser.add_option_group(runtimeSettings)
	
	# parse input files and/or metaphlan files from commandline
	try:
		options, x = parser.parse_args(argv)
	except:
		parser.error('[FATAL]: parsing failure.\n')
		exit(1)

	if not options.outdir:
		parser.error('[FATAL]: output directory is not supplied, you must specify output directory by -o/--outdir.\n')
		exit(1)
	
	if not options.config:
		parser.error('[FATAL]: configuration file is not supplied, you must specify the configuration file by -c/--config.\n')
		exit(1)
	
	script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
	
	if not options.ref_db:
		options.ref_db = script_path + '/db/ref_db'
		
	db_tarfile = options.ref_db + '.tar'
	db_tarindex = options.ref_db + '.tar.index'
	
	if not options.gsize_db:
		options.gsize_db = script_path + '/db/gsize.db'
		
	if not os.path.exists(db_tarfile):
		parser.error('[FATAL]: cannot locate species reference database tarball, %s' % db_tarfile)
	
	if not os.path.exists(db_tarindex):
		parser.error('[FATAL]: cannot locate species reference database index, %s' % db_tarindex)
	
	if not os.path.exists(options.gsize_db):
		parser.error('[FATAL]: cannot locate genome size database, %s' % options.gsize_db)
	
	if options.num_proc < 1:
		sys.stderr.write("[Warning]: illegal number of threads designated: %i, will use single thread if don\'t mind.\n" % options.num_proc)
		options.num_proc = 1
	
	if options.min_cov < 3:
		parser.error('[FATAL]: illegal minimum coverage supplied: %i; it has to be an int ranging in [3, +inf].\n' % options.min_cov)
	
	# parsing config file
	# parse config file
	samples = []	
	config_fh = open(options.config, 'r')
	sample = Sample()
	while 1:
		line = config_fh.readline().rstrip('\n')
		if not line: break
		if line[:2] == '//':
			if sample.is_empty(): continue
			else:
				samples.append(copy(sample))
				sample.clear()
				continue
		key, value = line.split(':')
		key = key.replace(' ', '', key.count(' '))
		key = key.replace('\t', '', key.count('\t'))
		value = value.replace(' ', '', value.count(' '))
		value = value.replace('\t', '', value.count('\t'))
		if key == 'sample':
			sample.name = value
		elif key == 'fq': 
			if os.path.exists(value): 
				sample.fq = os.path.abspath(value)
			else:
				sys.stderr.write('[FATAL]: Error parsing config file, encountered file non-existing: %s\n' % value)
				exit(1)
		elif key == 'fq1':
			if os.path.exists(value): 
				sample.fq1 = os.path.abspath(value)
			else:
				sys.stderr.write('[FATAL]: Error parsing config file, encountered file non-existing: %s\n' % value)
				exit(1)
		elif key == 'fq2':
			if os.path.exists(value): 
				sample.fq2 = os.path.abspath(value)
			else:
				sys.stderr.write('[FATAL]: Error parsing config file, encountered file non-existing: %s\n' % value)
				exit(1)
		elif key == 'metaphlan2':
			if os.path.exists(value): 
				sample.metaphlan2 = os.path.abspath(value)
			else:
				sys.stderr.write('[FATAL]: Error parsing config file, encountered file non-existing: %s\n' % value)
				exit(1)
		else:
			sys.stderr.write('[FATAL]: Error in keywords in config file: %s\n' % key)
			exit(1)
		
	if not sample.is_empty():
		samples.append(copy(sample))
		sample.clear()
	config_fh.close()
	
	# check if all required 3rd party programs are in place.
	if which(options.bowtie2) == None or which(options.bowtie2_build) == None:
		sys.stderr.write('[FATAL]: cannot locate Bowtie2 or Bowtie2-build.\n')
		exit(1)
	if which(options.samtools) == None:
		sys.stderr.write('[FATAL]: cannot locate samtools.\n')
		exit(1)
	if which(options.metaphlan2) == None:
		if need_metaphlan2(samples, options.outdir):
			sys.stderr.write('[FATAL]: cannot locate metaphlan2.\n')
			exit(1)
	
	# create output directory	
	if not os.path.exists(options.outdir): os.mkdir(options.outdir)
	
	# preserve project info into proj.log
	if not options.quiet: sys.stdout.write('\nProject info saved to: %s' % (options.outdir + '/proj.log'))
	proj_log = options.outdir + '/proj.log'
	proj_data = [options, samples, sys.argv]
	cPickle.dump(proj_data, open(proj_log, 'wb'))
	
	if not options.quiet: sys.stdout.write("\nNow constructing SNP flows...\n")
	construct(samples, db_tarfile, db_tarindex, options)
	if not options.quiet: sys.stdout.write("Done.\n")
	if not options.quiet: sys.stdout.write("Strain-GeMS construct finished!\n")

	## fitting the SNPs
	if not options.quiet: sys.stdout.write("\nNow fitting strains...\n")
	fit_strains(options, samples)
	if not options.quiet: sys.stdout.write("Done strain fitting!\n")
	
	## outputting results
	if not options.quiet: sys.stdout.write("\nNow outputting results...\n")
	output_results(samples, options)
	if not options.quiet: 
		sys.stdout.write("\nDone. Results are written to: %s.\n" % (options.outdir + '/results'))
		sys.stdout.write("Done. Strain-GeMS finished!\n")

if __name__ == '__main__':
	main()
