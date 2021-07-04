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
import os


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
  --dryrun                just show the commands we would have run, but don't
                          run them
  --treads=<N>            (T=) number of processing threads
                          (default is non-threaded)

Derive an estimate the proper size for the bloom filters in an SBT.

$$$$ more info here."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global ntcardCommand,howdesbtCommand
	global debug,isDryRun

	ntcardCommand   = "ntcard"
	howdesbtCommand = "howdesbt"

	# parse the command line

	fastqFilesFilename   = None
	kmerSize             = 20
	subsample            = None
	isDryRun             = False
	numThreads           = None
	candidateRatios      = None 
	overrideSizeEstimate = None
	debug                = []

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
		elif (arg == "--dryrun"):
			isDryRun = True
		elif (arg.startswith("--t=")) \
		  or (arg.startswith("--T=")) \
		  or (arg.startswith("--threads=")) \
		  or (arg.startswith("T=")):
			if (argVal == None):
				numThreads = None
			else:
				numThreads = int(argVal)
				if (numThreads <= 1): numThreads = None
		elif (arg.startswith("--candidate=")) or (arg.startswith("--candidates=")):  # (undocumented)
			if (candidateRatios == None): candidateRatios = []
			for candidateRatio in argVal.split(","):
				candidateRatio = float_or_fraction(candidateRatio)
				assert (0 < candidateRatio <= 2)
				candidateRatios += [candidateRatio]
		elif (arg.startswith("--sizeestimate=")):  # (undocumented)
			overrideSizeEstimate = int_with_unit(argVal)
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

	if (candidateRatios == None):
		candidateRatios = default_candidate_ratios()

	# read the fastq filenames

	filenamesF = open(fastqFilesFilename,"rt")
	fastqFilenames = read_filenames(filenamesF)
	filenamesF.close()

	assert (len(fastqFilenames) > 0), \
	       "\"%s\" contains no file names" % fastqFilenames

	for fastqFilename in fastqFilenames:
		assert (is_valid_fastq_name(fastqFilename)), \
		       "\"%s\" is not a valid fastq filename" % fastqFilename

	# get initial estimate the bloom filter size

	if (overrideSizeEstimate != None):
		bfSizeEstimate = overrideSizeEstimate
		print("using %s as simple bf size estimate" \
		    % "{:,}".format(bfSizeEstimate),
		      file=stderr)
	else:
		bfSizeEstimate = simple_bf_size_estimate(fastqFilenames,"temp.",kmerSize)
		print("simple bf size estimate is %s" \
		    % "{:,}".format(bfSizeEstimate),
		      file=stderr)

	# derive candidate sizes from the candidate expansion ratios (note that
	# many candidate ratios may be less than 1)

	bfSizeCandidates = [roundup64(r*bfSizeEstimate) for r in candidateRatios]
	print("bf size candidates are %s" \
	    % " ".join(["{:,}".format(bfSize) for bfSize in bfSizeCandidates]),
	      file=stderr)

	# generate the bloom filters for each candidate size

	if (numThreads == None):
		create_sbt_directories(bfSizeCandidates)
		for fastqFilename in fastqFilenames:
			howdesbt_make_bf(fastqFilename,kmerSize,bfSizeCandidates,subsample=subsample)

	else:
		assert (False), "threading is not implemented yet"

	# cluster the bloom filters for each candidate size

	if (numThreads == None):
		for numBits in bfSizeCandidates:
			howdesbt_cluster(fastqFilenames,numBits)

	else:
		assert (False), "threading is not implemented yet"

	# build the sequence bloom tree for each candidate size

	if (numThreads == None):
		for numBits in bfSizeCandidates:
			howdesbt_build(numBits)

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
	if (isDryRun) or ("ntcard" in debug):
		print("# running ntcard, output to \"%s\"\n%s" \
		    % (tempFilename," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		outcome = subprocess.run(command,capture_output=True)
		if (outcome.returncode != 0):
			dump_console_output(outcome)
			assert (False)


	if (isDryRun):
		return 1*1000*1000

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
		if (not isDryRun):
			outcome = subprocess.run(command,capture_output=True)

	# compute and return the estimate

	if ("ntcard" in debug) and (not isDryRun):
		print("ntcard abundance[F0] = %d" % abundance["F0"],file=stderr)
		print("ntcard abundance[f1] = %d" % abundance["1"],file=stderr)

	return abundance["F0"] - abundance["1"]


# create_sbt_directories--
#	Create directories to hold bloom filter files.

def create_sbt_directories(numBitsList):
	if (type(numBitsList) == int): numBitsList = [numBitsList]

	for numBits in numBitsList:
		command1 = ["rm","-rf",directory_name(numBits)]
		command2 = ["mkdir","-p",directory_name(numBits)]
		if (isDryRun) or ("howdesbt" in debug):
			print("# creating directory for numBits=%d\n%s\n%s" \
				% (numBits," ".join(command1)," ".join(command2)),
				  file=stderr)

		if (not isDryRun):
			outcome1 = subprocess.run(command1,capture_output=True)
			if (outcome1.returncode != 0):
				dump_console_output(outcome1)
				assert (False)

			outcome2 = subprocess.run(command2,capture_output=True)
			if (outcome2.returncode != 0):
				dump_console_output(outcome2)
				assert (False)


# howdesbt_make_bf--
#	Run HowDeSBT to create a bloom filter file (or several bloom filter files)
#	for a fastq file.

def howdesbt_make_bf(fastqFilename,kmerSize,numBitsList,subsample=None):
	if (type(numBitsList) == int): numBitsList = [numBitsList]
	fastqCoreName = fastq_core_name(fastqFilename)

	# ask howdesbt to build the bloom filter file(s)

	command = [howdesbtCommand,
	           "makebf",
	           "K=%d" % kmerSize,
	           "--min=2"]

	if (subsample == None):
		command += ["--bits=%d" % numBits for numBits in numBitsList]
	else:
		command += ["--modulus=%d" % numBits for numBits in numBitsList]
		command += ["--bits=%s" % subsample]

	command += [fastqFilename]

	if (subsample == None):
		command += ["--out=B={modulus}/%s.bf" % fastqCoreName]
	else:
		command += ["--out=B={bits}/%s.bf" % fastqCoreName]

	if (isDryRun) or ("howdesbt" in debug):
		print("# running howdesbt makebf for \"%s\"\n%s" \
		    % (fastqFilename," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		outcome = subprocess.run(command,capture_output=True)
		if (outcome.returncode != 0):
			dump_console_output(outcome)
			assert (False)


# howdesbt_cluster--
#	Run HowDeSBT to cluster bloom filter files as preparation for creating a
#	sequence bloom tree.

def howdesbt_cluster(fastqFilenames,numBits):
	workingDirectory  = directory_name(numBits)
	leafnamesFilename = "leafnames"
	treeFilename      = "union.sbt"

	# write the list of bloom filter files

	if (not isDryRun):
		f = open(workingDirectory+"/"+leafnamesFilename,"wt")
		for fastqCoreName in map(fastq_core_name,fastqFilenames):
			print("%s.bf" % fastqCoreName,file=f)
		f.close()

	# ask howdesbt to cluster the bloom filter files

	command = [howdesbtCommand,
	           "cluster",
	           "--list=%s" % leafnamesFilename,
	           "--tree=%s" % treeFilename,
	           "--nodename=node{number}",
	           "--keepallnodes"]

	if (isDryRun) or ("howdesbt" in debug):
		print("# running howdesbt cluster for %s\n%s\n%s" \
		    % (workingDirectory,"cd "+workingDirectory," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		outcome = subprocess.run(command,cwd=workingDirectory,capture_output=True)
		if (outcome.returncode != 0):
			dump_console_output(outcome)
			assert (False)

	if ("howdesbt" not in debug):
		command = ["rm",leafnamesFilename]
		if (not isDryRun): outcome = subprocess.run(command,capture_output=True)


# howdesbt_build--
#	Run HowDeSBT to build the sequence bloom tree.

def howdesbt_build(numBits):
	workingDirectory = directory_name(numBits)
	treeFilename     = "union.sbt"
	outTreeFilename  = "howde.sbt"

	# ask howdesbt to build the sequence bloom tree

	command = [howdesbtCommand,
	           "build",
	           "--HowDe",
	           "--tree=%s" % treeFilename,
	           "--outtree=%s" % outTreeFilename]

	if (isDryRun) or ("howdesbt" in debug):
		print("# running howdesbt build for %s\n%s\n%s" \
		    % (workingDirectory,"cd "+workingDirectory," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		outcome = subprocess.run(command,cwd=workingDirectory,capture_output=True)
		if (outcome.returncode != 0):
			dump_console_output(outcome)
			assert (False)


# directory_name--
#	Create a 'standardized' directory name for a particular bloom filter size
#	candidate.

def directory_name(numBits):
	return "B=%d" % numBits


# default_candidate_ratios--
#	Produce a list of candidate 'expansion' ratios that should be tried. These
#	are the ratio of bfSize/bfSizeEst, where bfSizeEst is a simple bloom filter
#	size estimate (e.g. from simple_bf_size_estimate), and bfSize is the size
#	of a candidate. Note that 'expansion' is misleading, as most of the ratios
#	are below 1.0.

def default_candidate_ratios():
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


# is_valid_fastq_name--
#	Determine whether a filename is one we recognize as a fastq file.

def is_valid_fastq_name(filename):
	return (filename.endswith(".fastq")) \
	    or (filename.endswith(".fq"))


# fastq_core_name--
#	Return the 'core' name of a fastq filename.

def fastq_core_name(filename):
	baseName = os.path.basename(filename)
	if (baseName.endswith(".fastq")): return baseName[:-len(".fastq")]
	if (baseName.endswith(".fq")):    return baseName[:-len(".fq")]
	raise ValueError


# dump_console_output--
#	Write the stdout and stderr output of the outcome of a call to 
#	subprocess.run that had capture_output=True.

def dump_console_output(outcome):
	for line in outcome.stdout.split(bytearray("\n",'utf-8')):
		print(line,file=stderr)
	for line in outcome.stderr.split(bytearray("\n",'utf-8')):
		print(line,file=stderr)


# roundup64--
#	Round a positive value up, to the next multiple of 64.

def roundup64(v):
	v = ceil(v)
	return ((v+63)//64)*64


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


# float_or_fraction--
#	Parse a value string, allowing a fraction, percentage, or floating point value

def float_or_fraction(s):
	if ("/" in s):
		(numer,denom) = s.split("/",1)
		return float(numer)/float(denom)
	elif (s.endswith("%")):
		return float(s[:1])/100
	else:
		return float(s)


if __name__ == "__main__": main()
