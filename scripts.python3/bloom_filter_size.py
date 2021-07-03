#!/usr/bin/env python3
"""
Derive an estimate the proper size for BFs in an SBT.

References:
  [1] Solomon, Brad, and Carl Kingsford. "Improved search of large
      transcriptomic sequencing databases using split sequence bloom trees."
      International Conference on Research in Computational Molecular Biology.
      Springer, Cham, 2017.
"""

from sys  import argv,stdin,stdout,stderr,exit
from math import sqrt,ceil
import subprocess


def usage(s=None):
	message = """

usage: bloom_filter_size <fastq_files_file> [options]
   <fastq_files_file>     file containing a list of fastq files to be
                          represented by the tree, one fastq file name per line
  --k=<N>                 (K=) kmer size (number of nucleotides in a kmer)
                          (default is 20)
  --subsample=<fraction>  fraction of bloom filter bits to be used for
                          evaluating a candidate size
                          (default is use the full candidate size)
  --treads=<N>            (T=) number of processing threads
                          (default is non-threaded)

Derive an estimate the proper size for the bloom filters in an SBT.

$$$$ more info here."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global ntcardCommand,howdesbtCommand
	global debug

	ntcardCommand   = "ntcard"
	howdesbtCommand = "howdesbt"

	# parse the command line

	fastqFilesFilename = None
	kmerSize           = 20
	subsample          = None
	numThreads         = None
	expansionRatios    = None   # $$$ add user options for this
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if   (arg.startswith("--k=")) \
		  or (arg.startswith("--K=")) \
		  or (arg.startswith("--kmer=")) \
		  or (arg.startswith("--kmersize=")) \
		  or (arg.startswith("K=")):
			kmerSize = int(argVal)
			assert (kmerSize > 0)
		elif (arg.startswith("--subsample=")):
			if (argVal == None):
				subsample = None
			else:
				subsample = float_or_fraction(argVal)
				if (subsample == 1):
					subsample = None
				else:
					assert (0 < subsample <= 1)
		elif (arg.startswith("--t=")) \
		  or (arg.startswith("--T=")) \
		  or (arg.startswith("--threads=")) \
		  or (arg.startswith("T=")):
			if (argVal == None):
				numThreads = None
			else:
				numThreads = int(argVal)
				if (numThreads <= 1): numThreads = None
		elif (arg.startswith("--ntcard=")):  # (undocumented)
			ntcardCommand = argVal
		elif (arg.startswith("--howdesbt=")):  # (undocumented)
			howdesbtCommand = argVal
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (fastqFilesFilename == None):
			fastqFilesFilename = arg
		else:
			usage("unrecognized option: %s" % arg)

	if (fastqFilesFilename == None):
		usage("you have to give me a fastq file names file")

	if (expansionRatios == None):
		expansionRatios = default_expansion_ratios()

	# read the fastq filenames

	filenamesF = open(fastqFilesFilename,"rt")
	fastqFilenames = read_filenames(filenamesF)
	filenamesF.close()

	assert (len(fastqFilenames) > 0), \
	       "\"%s\" contains no file names" % fastqFilenames

	# get initial estimate the bloom filter size

	bfSizeEstimate = simple_bf_size_estimate(fastqFilenames,"temp.",kmerSize)
	print("simple bf size estimate is %s" \
	    % "{:,}".format(bfSizeEstimate),
	      file=stderr)

	# derive candidate sizes from the expansion ratios (note that many
	# expansion ratios may be less than 1)

	bfSizeCandidates = [roundup64(r*bfSizeEstimate) for r in expansionRatios]
	print("bf size candidates are %s" \
	    % " ".join(["{:,}".format(bfSize) for bfSize in bfSizeCandidates]),
	      file=stderr)

	# generate the candidate bloom filters

	if (numThreads == None):
		for fastqFilename in fastqFilenames:
			howdesbt_make_bf(fastqFilename,kmerSize,bfSizeCandidates,subsample=subsample)

	else:
		assert (False), "threading is not implemented yet"


# simple_bf_size_estimate--
#	Derive a simple estimate for the size of the bloom filters. This follows
#	the recommendation in Solomon2017 [1], which is to use the number of
#	distinct kmers among all the experiments. Here we assume that single kmers
#	will be excluded from the tree, so we want F0-f1, where F0 is the number of
#	distinct kmers and f1 is the number of singletons.
#
# Output from ntcard will look like this:
#	F1  7806530
#	F0  2114946
#	1   1360743     <--- this is f1
#	2   74131
#	3   73667
#    ...
# Earlier versions of ntcard output a lowercase f prefix on all the lines that
# now have no prefix (e.g. "1 1360743" used to be "f1 1360743"), so we tolerate
# that format.

def simple_bf_size_estimate(fastqFilenames,tempFilenamePrefix,kmerSize):
	# ask ntcard to count kmer abundances

	tempFilename = "%s_k%d.hist" % (tempFilenamePrefix,kmerSize)
	command = [ntcardCommand,
	           "--kmer=%d" % kmerSize,
	           "--pref=%s" % tempFilenamePrefix] \
	        + fastqFilenames
	if ("ntcard" in debug):
		print("running ntcard, output to \"%s\"\n%s" \
		    % (tempFilename," ".join(command)),
		      file=stderr)
	outcome = subprocess.run(command,capture_output=True)
	assert (outcome.returncode == 0)

	# parse ntcard's output

	abundance = {}

	f = open(tempFilename,"rt")

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		fields = line.split()

		try:
			if (len(fields) != 2): raise ValueError
			key   = fields[0]
			count = int(fields[1][1:]) if (fields[1].startswith("f")) else int(fields[1])
			abundance[key] = count
		except ValueError:
			assert (len(fields) == 2), \
			       "unexpected ntcard output at line %d of \"%s\":\n%s" \
			     % (lineNumber,tempFilename,line)

	f.close()

	assert ("F0" in abundance) and ("1" in abundance), \
	       "ntcard output lacks a value for F0 and/or 1, (see \"%s\")" \
	     % tempFilename

	if ("ntcard" not in debug):
		command = ["rm",tempFilename]
		outcome = subprocess.run(command,capture_output=True)

	# compute and return the estimate

	if ("ntcard" in debug):
		print("ntcard abundance[F0] = %d" % abundance["F0"],file=stderr)
		print("ntcard abundance[f1] = %d" % abundance["1"],file=stderr)

	return abundance["F0"] - abundance["1"]


# howdesbt_make_bf--
#	Run HowDeSBT to create a bloom filter file (or several bloom filter files)
#	for a fastq file.

def howdesbt_make_bf(fastqFilename,kmerSize,numBits,subsample=None):
	# ask howdesbt to build the bloom filter file(s)

#$$$ change this
	if (type(numBits) != int): numBits = numBits[0]

	command = [howdesbtCommand,
	           "makebf",
	           "K=%d" % kmerSize,
	           "--min=2",
	           "--bits=%d" % numBits,
	           fastqFilename]
	if ("howdesbt" in debug):
		print("running howdesbt for \"%s\"\n%s" \
		    % (fastqFilename," ".join(command)),
		      file=stderr)
	outcome = subprocess.run(command,capture_output=True)
	assert (outcome.returncode == 0)


# default_expansion_ratios--
#	Produce a list of candidate 'expansion' ratios that should be tried. These
#	are the ratio of bfSize/bfSizeEst, where bfSizeEst is a simple bloom filter
#	size estimate (e.g. from simple_bf_size_estimate), and bfSize is the size
#	of a candidate. Note that 'expansion' is misleading, as most of the ratios
#	are below 1.0.

def default_expansion_ratios():
	r = sqrt(5/4)
	return [r**n for n in range(-10,3)]


# read_filenames--
#	Read filenames from a text file, one per line, stripping any whitespace.

def read_filenames(f):
	filenames = []
	for line in f:
		filename = line.strip()
		if (filename == ""): continue
		filenames += [filename]
	return filenames


# roundup64--
#	Round a positive value up, to the next multiple of 64.

def roundup64(v):
	v = ceil(v)
	return ((v+63)//64)*64


# float_or_fraction--
#	Parse a value string, allowing a fraction, percentage, or floating point value

def float_or_fraction(s):
	if ("/" in s):
		(numer,denom) = s.split("/",1)
		return float(numer)/float(denom)
	elif (s.endswith("%")):
		float(s[:1])/100
	else:
		return float(s)


if __name__ == "__main__": main()
