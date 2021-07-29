#!/usr/bin/env python3
"""
Derive an estimate the proper size for BFs in an SBT.

References:
  [1] Harris, Robert S., and Paul Medvedev. "Improved representation of
      sequence bloom trees." Bioinformatics 36.3 (2020): 721-727.
  [2] Solomon, Brad, and Carl Kingsford. "Improved search of large
      transcriptomic sequencing databases using split sequence bloom trees."
      International Conference on Research in Computational Molecular Biology.
      Springer, Cham, 2017.
"""

from sys  import argv,stdin,stdout,stderr,exit,version_info
from math import sqrt,log,ceil
from multiprocessing.pool import ThreadPool as Pool
import subprocess
import os
from stat import S_IEXEC

standardStep = sqrt(5/4)


def usage(s=None):
	message = """

usage: bloom_filter_size <input_files_file> [options]
   <input_files_file>     file containing a list of fastq (or other) files to
                          be represented by the tree, one fastq file name per
                          line
  --k=<N>                 (K=) kmer size (number of nucleotides in a kmer)
                          (default is 20)
  --stable=<percent>      stabilization threshold (see note below)
                          (default is 5%)
  --subsample=<fraction>  fraction of bloom filter bits to be used for
                          evaluating a candidate size
                          (default is use the 1/40th of the candidate size)
  --subsample=none        use the full candidate size for evaluating a
                          candidate size
  --cluster=<fraction>    fraction of bloom filter bits to be used for
                          clustering
                          (default is use the full candidate or subsample size)
  --dryrun                just show the commands we would have run, but don't
                          run them
  --threads=<N>           (T=) number of processing threads
                          (default is non-threaded)

Derive an estimate the proper size for the bloom filters in an SBT.

The <input_files_file> contains a list of files, each to represented by a leaf
int the SBT. These are usually fast files, but jellyfish kmer-count files are
also allowed. Valid extensions are .fastq, .fq, .fastq.gz, .fq.gz, .jellyfish,
.jf, .jellyfish.gz, or .jf.gz.

The 'stabilization threshold' is the relative increase in tree size (as bytes
on disk), below which we consider the size to be stable. The BF size we report
(as BF vector size in bits) is the lowest candidate which is stable. The
increase in tree size is measured relative to an increase in BF size of
sqrt(5/4)."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global ntcardCommand,howdesbtCommand
	global debug,isDryRun

	ntcardCommand   = "ntcard"
	howdesbtCommand = "howdesbt"

	# parse the command line

	fastqFilesFilename     = None
	kmerSize               = 20
	stabilizationThreshold = .05
	subsampleFraction      = 1./40
	clusterFraction        = None
	isDryRun               = False
	numThreads             = None
	candidateRatios        = None 
	candidateExponents     = None
	overrideSizeEstimate   = None
	debug                  = []

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
		elif (arg.startswith("--stable=")):
			stabilizationThreshold = float_or_fraction(argVal)
			assert (0 < stabilizationThreshold < 1)
		elif (arg.startswith("--subsample=")):
			if (argVal == "none"):
				subsampleFraction = None
			else:
				subsampleFraction = float_or_fraction(argVal)
				if (subsampleFraction == 1):
					subsampleFraction = None
				else:
					assert (0 < subsampleFraction <= 1)
		elif (arg == "--nosubsample"):
			subsampleFraction = None
		elif (arg.startswith("--cluster=")):
			if (argVal == None):
				clusterFraction = None
			else:
				clusterFraction = float_or_fraction(argVal)
				if (clusterFraction == 1):
					clusterFraction = None
				else:
					assert (0 < clusterFraction <= 1)
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
		elif (arg.startswith("--exponents=")):  # (undocumented)
			(start,end) = argVal.split(",",2)
			candidateExponents = (int(start),int(end))
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
		candidateRatios = default_candidate_ratios(exponents=candidateExponents)

	# read the fastq filenames

	filenamesF = open(fastqFilesFilename,"rt")
	inputFilenames = read_filenames(filenamesF)
	filenamesF.close()

	assert (len(inputFilenames) > 0), \
	       "\"%s\" contains no file names" % inputFilenames

	for inputFilename in inputFilenames:
		assert (input_file_type(inputFilename) != None), \
		       "\"%s\" is not a valid fastq/input filename" % inputFilename
		if (overrideSizeEstimate == None):
			assert (input_file_type(inputFilename) in ["fastq"]), \
			       "\"%s\" is not a valid fastq/input filename for ntcard" % inputFilename

	# get initial estimate the bloom filter size

	if (overrideSizeEstimate != None):
		bfSizeEstimate = overrideSizeEstimate
		print("using %s as simple bf size estimate" \
		    % "{:,}".format(bfSizeEstimate),
		      file=stderr)
	else:
		bfSizeEstimate = simple_bf_size_estimate(inputFilenames,"temp.",kmerSize)
		print("simple bf size estimate is %s" \
		    % "{:,}".format(bfSizeEstimate),
		      file=stderr)

	# derive candidate sizes from the candidate expansion ratios (note that
	# many candidate ratios may be less than 1)

	# was bfSizeCandidates = [roundup64(r*bfSizeEstimate) for r in candidateRatios]
	bfSizeCandidates = [roundup(r*bfSizeEstimate) for r in candidateRatios]
	bfSizeCandidates.sort()
	print("bf size candidates are %s" \
	    % " ".join(["{:,}".format(bfSize) for bfSize in bfSizeCandidates]),
	      file=stderr)

	# generate the uncompressed bloom filters for each candidate size

	create_sbt_directories(bfSizeCandidates)

	jobs = []
	for (num,inputFilename) in enumerate(inputFilenames):
		job = SbtMakeBfJob()
		job.inputFilename     = inputFilename
		job.kmerSize          = kmerSize
		job.numBitsList       = bfSizeCandidates
		job.subsampleFraction = subsampleFraction
		job.comment           = "(#%d of %d)" % (num+1,len(inputFilenames))
		jobs += [job]

	if (numThreads == None):
		for _ in map(howdesbt_make_bf,jobs):
			pass # $$$ fix this
	else:
		pool = Pool(numThreads)
		for _ in pool.imap_unordered(howdesbt_make_bf,jobs):
			pass # $$$ fix this

	# cluster the bloom filters for each candidate size

	jobs = []
	for (num,numBits) in enumerate(bfSizeCandidates):
		job = SbtClusterJob()
		job.inputFilenames  = inputFilenames
		job.numBits         = numBits
		job.clusterFraction = clusterFraction
		job.comment         = "(#%d of %d)" % (num+1,len(bfSizeCandidates))
		jobs += [job]

	if (numThreads == None):
		for _ in map(howdesbt_cluster,jobs):
			pass # $$$ fix this
	else:
		pool = Pool(numThreads)
		for _ in pool.imap_unordered(howdesbt_cluster,jobs):
			pass # $$$ fix this

	# build the sequence bloom tree for each candidate size

	jobs = []
	for (num,numBits) in enumerate(bfSizeCandidates):
		job = SbtBuildJob()
		job.numBits        = numBits
		job.inputFilenames = inputFilenames
		job.comment        = "(#%d of %d)" % (num+1,len(inputFilenames))
		jobs += [job]

	if (numThreads == None):
		for _ in map(howdesbt_build,jobs):
			pass # $$$ fix this
	else:
		pool = Pool(numThreads)
		for _ in pool.imap_unordered(howdesbt_build,jobs):
			pass # $$$ fix this

	# record the total size-on-disk of each sequence bloom tree

	jobs = []
	for (num,numBits) in enumerate(bfSizeCandidates):
		job = SizeOnDiskJob()
		job.numBits = numBits
		jobs += [job]

	numBitsToSizeOnDisk = {}
	if (numThreads == None):
		for (numBits,sizeOnDisk) in map(size_on_disk,jobs):
			numBitsToSizeOnDisk[numBits] = sizeOnDisk
	else:
		pool = Pool(numThreads)
		jobOrder = []
		for (numBits,sizeOnDisk) in pool.imap_unordered(size_on_disk,jobs):
			numBitsToSizeOnDisk[numBits] = sizeOnDisk
			jobOrder += [numBits]

		if ("pool" in debug):
			print("job order: %s" % ",".join(["%d"%numBits for numBits in jobOrder]),
			      file=stderr)

	# choose the threshold

	meetsThreshold        = {}
	numBitsToNormStep     = {}
	numBitsToIncrease     = {}
	numBitsToNormIncrease = {}

	for (ix,numBits) in enumerate(bfSizeCandidates):
		if (ix == 0): continue
		prevNumBits = bfSizeCandidates[ix-1]
		assert (numBits > prevNumBits)

		prevSizeOnDisk = numBitsToSizeOnDisk[prevNumBits]
		sizeOnDisk     = numBitsToSizeOnDisk[numBits]

		if (sizeOnDisk < prevSizeOnDisk):
			meetsThreshold[numBits] = True
			continue

		step               = numBits / prevNumBits
		normalizedStep     = log(step) / log(standardStep)
		increase           = (sizeOnDisk / prevSizeOnDisk) - 1
		normalizedIncrease = ((1+increase)**(1/normalizedStep)) - 1

		meetsThreshold[numBits]        = (normalizedIncrease <= stabilizationThreshold)
		numBitsToNormStep[numBits]     = normalizedStep
		numBitsToIncrease[numBits]     = increase
		numBitsToNormIncrease[numBits] = normalizedIncrease

	# report results

	print("#%s\t%s\t%s\t%s\t%s\t%s\t%s" \
	    % ("bfBits","step","sizeOnDisk","predictedSize",
	       "increase%","normalized%","stable"))

	firstStableNumBits = None

	for numBits in bfSizeCandidates:
		sizeOnDisk = numBitsToSizeOnDisk[numBits]

		if (firstStableNumBits == None):
			if (numBits in meetsThreshold) and (meetsThreshold[numBits]):
				firstStableNumBits = numBits

		print("%d\t%s\t%d\t%d\t%s\t%s\t%s" \
			% (numBits,
			   "NA" if (numBits not in numBitsToNormStep) \
					else ("%0.2f" % numBitsToNormStep[numBits]),
			   sizeOnDisk,
			   int(ceil(sizeOnDisk/subsampleFraction)),
			   "NA" if (numBits not in numBitsToIncrease) \
					else ("%0.2f%%" % (100*numBitsToIncrease[numBits])),
			   "NA" if (numBits not in numBitsToNormIncrease) \
					else ("%0.2f%%" % (100*numBitsToNormIncrease[numBits])),
			   "NA" if (numBits not in meetsThreshold) \
					else ["no","yes"][meetsThreshold[numBits]]))

	if (firstStableNumBits != None):
		print("proper size is estimated to be %d" % firstStableNumBits)
	else:
		print("no size met the stability threshold, max tried was %d" % max(bfSizeCandidates))

	# $$$ should we cleanup, discarding all the candidate bf directories?


# simple_bf_size_estimate--
#	Derive a simple estimate for the size of the bloom filters. This follows
#	the recommendation in Solomon2017 [2], which is to use the number of
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

def simple_bf_size_estimate(inputFilenames,tempFilenamePrefix,kmerSize):
	# ask ntcard to count kmer abundances

	tempFilename = "%s_k%d.hist" % (tempFilenamePrefix,kmerSize)
	command = [ntcardCommand,
	           "--kmer=%d" % kmerSize,
	           "--pref=%s" % tempFilenamePrefix] \
	        + inputFilenames
	if (isDryRun) or ("ntcard" in debug):
		print("# running ntcard, output to \"%s\"\n%s" \
		    % (tempFilename," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		run_shell_command(command)

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

		#if ("ntcard" in debug):
		#	print("line %d: \"%s\"" % (lineNumber,fields),file=stderr)

		try:
			if (len(fields) != 2): raise ValueError
			if   (fields[0].startswith("F")): key = fields[0]
			elif (fields[0].startswith("f")): key = int(fields[0][1:])
			else:                             key = int(fields[0])
			count = int(fields[1])
			abundance[key] = count
		except ValueError:
			assert (len(fields) == 2), \
			       "unexpected ntcard output at line %d of \"%s\":\n%s" \
			     % (lineNumber,tempFilename,line)

	f.close()

	assert ("F0" in abundance) and (1 in abundance), \
	       "ntcard output lacks a value for F0 and/or 1, (see \"%s\")" \
	     % tempFilename

	if ("ntcard" not in debug):
		command = ["rm",tempFilename]
		if (not isDryRun):
			run_shell_command(command)

	# compute and return the estimate

	if ("ntcard" in debug) and (not isDryRun):
		print("ntcard abundance[F0] = %d" % abundance["F0"],file=stderr)
		print("ntcard abundance[f1] = %d" % abundance[1],file=stderr)

	return abundance["F0"] - abundance[1]


# create_sbt_directories--
#	Create directories to hold bloom filter files.

def create_sbt_directories(numBitsList):
	if (type(numBitsList) == int): numBitsList = [numBitsList]

	for numBits in numBitsList:
		command1 = ["rm","-rf",candidate_directory_name(numBits)]
		command2 = ["mkdir","-p",candidate_directory_name(numBits)]
		if (isDryRun) or ("howdesbt" in debug):
			print("# creating directory for numBits=%d\n%s\n%s" \
				% (numBits," ".join(command1)," ".join(command2)),
				  file=stderr)

		if (not isDryRun):
			run_shell_command(command1)
			run_shell_command(command2)


# howdesbt_make_bf--
#	Run HowDeSBT to create a bloom filter file (or several bloom filter files)
#	for a fastq file.
#
# Implementation notes:
#	I was unable to implement a deadlock free version of run_shell_command for
#	commands with pipelines. Thus I have opted to write my command to a
#	temporary shell script and run that.

class SbtMakeBfJob: pass

def howdesbt_make_bf(job):
	if (type(job.numBitsList) == int): numBitsList = [job.numBitsList]
	else:                              numBitsList = list(job.numBitsList)

	fileType = input_file_type(job.inputFilename)
	inputCoreName = input_core_name(job.inputFilename)
	candidateDirectory = candidate_directory_template(job.subsampleFraction)
	comment = "" if (job.comment == None) else (" "+job.comment)

	minAbundance = 2

	# ask howdesbt to build the bloom filter file(s)

	if (job.subsampleFraction == None):
		bitsArgs = ["--bits=%d" % numBits for numBits in numBitsList]
	else:
		bitsArgs = ["--modulus=%d" % numBits for numBits in numBitsList] \
		         + ["--bits=%s%%" % (100*job.subsampleFraction)]

	if (fileType == "gzipped fastq"):
		command = ["gzip","-dc",job.inputFilename,"|",
		           howdesbtCommand,"makebf",
		           "K=%d" % job.kmerSize,
		           "--min=%d" % minAbundance] \
		        + bitsArgs \
		        + ["/dev/stdin",
		           "--out=%s/%s.bf" % (candidateDirectory,inputCoreName)]
	elif (fileType == "jellyfish"):
		command = ["jellyfish","dump","--column","--lower-count=%d"%minAbundance,job.inputFilename,"|",
		           howdesbtCommand,"makebf",
		           "--kmersin","K=%d" % job.kmerSize] \
		        + bitsArgs \
		        + ["/dev/stdin",
		           "--out=%s/%s.bf" % (candidateDirectory,inputCoreName)]
	elif (fileType == "gzipped jellyfish"):
		command = ["gzip","-dc",job.inputFilename,"|",
		           "jellyfish","dump","--column","--lower-count=%d"%minAbundance,"/dev/stdin","|",
		           howdesbtCommand,"makebf",
		           "--kmersin","K=%d" % job.kmerSize] \
		        + bitsArgs \
		        + ["/dev/stdin",
		           "--out=%s/%s.bf" % (candidateDirectory,inputCoreName)]
	else:
		command = [howdesbtCommand,"makebf",
		           "K=%d" % job.kmerSize,
		           "--min=%d" % minAbundance] \
		        + bitsArgs \
		        + [job.inputFilename,
		           "--out=%s/%s.bf" % (candidateDirectory,inputCoreName)]

	# if the command is a pipeline, write it as a temporary shell script

	scriptFilename = None
	withShell = False

	commandIsPipeline = ("|" in command)
	if (commandIsPipeline):
		scriptDirectory = candidate_directory_name(numBitsList[0])
		scriptFilename = "%s/%s.sh" % (scriptDirectory,inputCoreName)
		f = open(scriptFilename,"wt")
		print(" ".join(command),file=f)
		f.close()
		os.chmod(scriptFilename,os.stat(scriptFilename).st_mode|S_IEXEC) # (make it executable)

		if (isDryRun) or ("howdesbt" in debug):
			print("# creating temporary shell script \"%s\"\n%s" \
			    % (scriptFilename," ".join(command)),
			      file=stderr)

		if ("commands" in debug):
			print("[%s] %s" % (scriptFilename," ".join(command)),file=stderr)

		withShell = True
		command = [scriptFilename]

	# run the command

	if (isDryRun) or ("howdesbt" in debug):
		print("#%s running howdesbt makebf for \"%s\"\n%s" \
		    % (comment,job.inputFilename," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		run_shell_command(command,withShell=withShell)

	if (scriptFilename != None):
		command = ["rm",scriptFilename]
		if (not isDryRun):
			run_shell_command(command)

	return None


# howdesbt_cluster--
#	Run HowDeSBT to cluster bloom filter files as preparation for creating a
#	sequence bloom tree.

class SbtClusterJob: pass

def howdesbt_cluster(job):
	workingDirectory  = candidate_directory_name(job.numBits)
	leafnamesFilename = "leafnames"
	treeFilename      = "union.sbt"
	comment = "" if (job.comment == None) else (" "+job.comment)

	# write the list of bloom filter files

	if (not isDryRun):
		f = open(workingDirectory+"/"+leafnamesFilename,"wt")
		for inputCoreName in map(input_core_name,job.inputFilenames):
			print("%s.bf" % inputCoreName,file=f)
		f.close()

	# ask howdesbt to cluster the bloom filter files

	command =  [howdesbtCommand,
	            "cluster",
	            "--list=%s" % leafnamesFilename]
	if (job.clusterFraction != None):
		# was command += ["--bits=%d" % roundup64(job.clusterFraction*job.numBits)]
		command += ["--bits=%d" % roundup(job.clusterFraction*job.numBits)]
	command += ["--tree=%s" % treeFilename,
	            "--nodename=node{number}",
	            "--keepallnodes"]

	if (isDryRun) or ("howdesbt" in debug):
		print("#%s running howdesbt cluster for %s\n%s\n%s" \
		    % (comment,workingDirectory,"cd "+workingDirectory," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		run_shell_command(command,cwd=workingDirectory)

	if ("howdesbt" not in debug):
		command = ["rm",leafnamesFilename]
		if (not isDryRun):
			run_shell_command(command,cwd=workingDirectory)


# howdesbt_build--
#	Run HowDeSBT to build a sequence bloom tree.

class SbtBuildJob: pass

def howdesbt_build(job):
	workingDirectory = candidate_directory_name(job.numBits)
	treeFilename     = "union.sbt"
	outTreeFilename  = "howde.sbt"
	comment = "" if (job.comment == None) else (" "+job.comment)

	# ask howdesbt to build the sequence bloom tree

	command = [howdesbtCommand,
	           "build",
	           "--HowDe",
	           "--tree=%s" % treeFilename,
	           "--outtree=%s" % outTreeFilename]

	if (isDryRun) or ("howdesbt" in debug):
		print("#%s running howdesbt build for %s\n%s\n%s" \
		    % (comment,workingDirectory,"cd "+workingDirectory," ".join(command)),
		      file=stderr)
	if (not isDryRun):
		run_shell_command(command,cwd=workingDirectory)

	# delete the input tree file and the uncompressed bloom filters

	if ("howdesbt" not in debug):
		command = ["rm","%s/%s" % (workingDirectory,treeFilename)]
		if (not isDryRun):
			run_shell_command(command)
		for inputFilename in job.inputFilenames:
			inputCoreName = input_core_name(inputFilename)
			command = ["rm","%s/%s.bf" % (workingDirectory,inputCoreName)]
			if (not isDryRun):
				run_shell_command(command)


# size_on_disk--
#	Collect the size-on-disk of all the files in a sequence bloom tree.
#
# For information about using os.stat() to compute file size, see
#	https://stackoverflow.com/questions/2104080/how-can-i-check-file-size-in-python

class SizeOnDiskJob: pass

def size_on_disk(job):
	workingDirectory = candidate_directory_name(job.numBits)
	treeFilename     = "howde.sbt"

	if (isDryRun):
		print("# getting sizes for %s" % workingDirectory+"/"+treeFilename,file=stderr)
		return 0

	f = open(workingDirectory+"/"+treeFilename,"rt")
	bfFilenames = read_sbt_filenames(f)
	f.close()

	totalSizeOnDisk = 0
	for bfFilename in bfFilenames:
		path = workingDirectory + "/" + bfFilename
		sizeOnDisk = os.stat(path).st_size
		totalSizeOnDisk += sizeOnDisk
		if ("sizeondisk" in debug):
			print("size: %s %d" % (path,sizeOnDisk),file=stderr)

	return (job.numBits,totalSizeOnDisk)


# candidate_directory_name, candidate_directory_template--
#	Create a 'standardized' directory template, or directory name for a
#	particular bloom filter size candidate.

def candidate_directory_template(subsampleFraction):
	if (subsampleFraction == None):
		return "temp.B={bits}"
	else:
		return "temp.B={modulus}"

def candidate_directory_name(numBits):
	return "temp.B=%d" % numBits


# default_candidate_ratios--
#	Produce a list of candidate 'expansion' ratios that should be tried. These
#	are the ratio of bfSize/bfSizeEst, where bfSizeEst is a simple bloom filter
#	size estimate (e.g. from simple_bf_size_estimate), and bfSize is the size
#	of a candidate. Note that 'expansion' is misleading, as most of the ratios
#	are below 1.0.

def default_candidate_ratios(exponents=None):
	if (exponents == None): exponents = (-10,3)
	(start,end) = exponents
	r = standardStep
	return [r**n for n in range(start,end)]


# run_shell_command--
#	Run a command in the shell, capturing its console output.
#
# Implementation notes:
#	(1)	subprocess.run(...capture_output...) requires python 3.7 or newer
#	(2)	As best I can tell, using subprocess.run(shell=True) to support the
#	    pipeline case suffers from deadlock when one command outputs too much
#	    to stdout. 
#	(3) Though the pipeline case is briefly discussed in the python docs at
#		  docs.python.org/3/library/subprocess.html#replacing-shell-pipeline
#	    and a three-command example appears at
#	      stackoverflow.com/questions/295459/how-do-i-use-subprocess-popen-to-connect-multiple-processes-by-pipes
#	    I was unable to implement a deadlock free version.

class Outcome: pass

def run_shell_command(command,withShell=False,cwd=None):
	assert (type(command) == list)
	commandIsPipeline = ("|" in command)

	if ("commands" in debug):
		print("[command] "+" ".join(command),file=stderr)

	assert (not commandIsPipeline)

	if (version_info >= (3,7)): # (python3.7 or newer)
		outcome = subprocess.run(command,shell=withShell,cwd=cwd,capture_output=True)
	else:
		outcome = subprocess.run(command,shell=withShell,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

	if (outcome.returncode != 0):
		dump_console_output(command,outcome)
		assert (False)

	del outcome.stdout
	del outcome.stderr

	return outcome


# read_filenames--
#	Read filenames from a text file, one per line, stripping any whitespace.

def read_filenames(f):
	filenames = []
	for line in f:
		filename = line.strip()
		if (filename == ""): continue
		filenames += [filename]
	return filenames


# read_sbt_filenames--
#	Read bloom filter filenames from a file in sbt format.
#
# The sbt tree hierarchy file looks something like this (below). For the
# purpose of this function, if consists of a filename and an optional prefix of
# asterisks (which is ignored).
#	node1.bf
#	*node2.bf
#	**node4.bf
#	***node8.bf
#	****EXPERIMENT18.bf
#	****EXPERIMENT7.bf
#	***node9.bf
#	****EXPERIMENT12.bf
#	 ...

def read_sbt_filenames(f):
	filenames = []
	for line in f:
		line = line.strip()
		indent = 0
		while (line[indent] == "*"):
			indent += 1
		filename = line[indent:]
		filenames += [filename]
	return filenames


# input_file_type--
#	Determine what type of input file a filename represents, and whether it is
#	one that we support.
#
# Note that we reject filenames containing any whitespace, so that we don't
# have to worry about whether the commands we use support whitespace.

def input_file_type(filename):
	if (filename != "".join(filename.split())):
		return None
	if (filename.endswith(".fastq")) or (filename.endswith(".fq")):
		return "fastq"
	if (filename.endswith(".fastq.gz")) or (filename.endswith(".fq.gz")):
		return "gzipped fastq"
	if (filename.endswith(".jellyfish")) or (filename.endswith(".jf")):
		return "jellyfish"
	if (filename.endswith(".jellyfish.gz")) or (filename.endswith(".jf.gz")):
		return "gzipped jellyfish"
	return None


# input_core_name--
#	Return the 'core' name of a fastq or jellyfish filename.

def input_core_name(filename):
	baseName = os.path.basename(filename)
	if (baseName.endswith(".fastq")):        return baseName[:-len(".fastq")]
	if (baseName.endswith(".fq")):           return baseName[:-len(".fq")]
	if (baseName.endswith(".fastq.gz")):     return baseName[:-len(".fastq.gz")]
	if (baseName.endswith(".fq.gz")):        return baseName[:-len(".fq.gz")]
	if (baseName.endswith(".jellyfish")):    return baseName[:-len(".jellyfish")]
	if (baseName.endswith(".jf")):           return baseName[:-len(".jf")]
	if (baseName.endswith(".jellyfish.gz")): return baseName[:-len(".jellyfish.gz")]
	if (baseName.endswith(".jf.gz")):        return baseName[:-len(".jf.gz")]
	raise ValueError


# dump_console_output--
#	Write the stdout and stderr output of the outcome of a call to 
#	subprocess.run that had capture_output=True.

# $$$ make this print strings instread of byte arrays

def dump_console_output(command,outcome):
	isFirst = True
	for line in outcome.stdout.decode().split("\n"):
		if (isFirst):
			print(" ".join(command),file=stderr)
			isFirst = False
		print("[stdout] "+line,file=stderr)
	for line in outcome.stderr.decode().split("\n"):
		if (isFirst):
			print(" ".join(command),file=stderr)
			isFirst = False
		print("[stderr] "+line,file=stderr)


# roundup--
#	Round a positive value up, to the next integer

def roundup(v):
	return int(ceil(v))


# roundup64--
#	Round a positive value up, to the next multiple of 64.

def roundup64(v):
	v = int(ceil(v))
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
		return float(s[:-1])/100
	else:
		return float(s)


if __name__ == "__main__": main()
