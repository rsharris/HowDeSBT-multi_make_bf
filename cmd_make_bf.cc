// cmd_make_bf.cc-- convert a sequence file to a bloom filter

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>

#include "utilities.h"
#include "hash.h"
#include "jelly_kmers.h"
#include "bloom_filter.h"
#include "bloom_filter_file.h"

#include "support.h"
#include "commands.h"
#include "cmd_make_bf.h"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
#define u32 std::uint32_t
#define u64 std::uint64_t

void MakeBFCommand::short_description
   (std::ostream& s)
	{
	s << commandName << "-- convert a sequence file to a bloom filter" << endl;
	}

void MakeBFCommand::usage
   (std::ostream& s,
	const string& message)
	{
	if (!message.empty())
		{
		s << message << endl;
		s << endl;
		}

	short_description(s);
	s << "usage: " << commandName << " <filename> [<filename>..] [options]" << endl;
	//    123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	s << "  <filename>         (cumulative) a sequence file, e.g. fasta, fastq, or kmers" << endl;
	s << "                     (one bloom filter is created, for the union of the" << endl;
	s << "                     sequence files)" << endl;
	s << "  --kmersin          input files are kmers" << endl;
	s << "                     (by default input files are expected to be fasta or fastq)" << endl;
	s << "  --out=<filename>   name for bloom filter file" << endl;
	s << "                     (by default this is derived from first sequence filename)" << endl;
	s << "  --list=<filename>  file containing a list of bloom filters to create; this is" << endl;
	s << "                     used in place of the <filename>s on the command line; the" << endl;
	s << "                     file format is described below" << endl;
	s << "  --asper=<filename> name of an existing bloom filter file to extract settings" << endl;
	s << "                     from; that file's --k, --hashes, --seed, --modulus," << endl;
	s << "                     --bits and compression type will be used if they are not" << endl;
	s << "                     otherwise specified on the command line" << endl;
	s << "  --k=<N>            kmer size (number of nucleotides in a kmer)" << endl;
	s << "                     (default is " << defaultKmerSize << ")" << endl;
	s << "  --min=<N>          kmers occuring fewer than N times are left out of the" << endl;
	s << "                     bloom filter; this does not apply when --kmersin is used" << endl;
	s << "                     (default is " << defaultMinAbundance << ")" << endl;
	s << "  --threads=<N>      number of threads to use during kmerization" << endl;
	s << "                     (default is " << defaultNumThreads << ")" << endl;
	s << "  --hashes=<N>       how many hash functions to use for the filter" << endl;
	s << "                     (default is " << defaultNumHashes << ")" << endl;
	s << "  --seed=<number>    the hash function's 56-bit seed" << endl;
	s << "  --seed=<number>,<number>  both the hash function seeds; the second seed is" << endl;
	s << "                     only used if more than one hash function is being used" << endl;
	s << "                     (by default the second seed is the first seed plus 1)" << endl;
	s << "  --modulus=<M>      (cumulative) set the hash modulus, if larger than the" << endl;
	s << "                     number of bits" << endl;
	s << "                     (by default this is the same as the number of bits)" << endl;
	s << "  --bits=<N>         (cumulative) number of bits in the bloom filter" << endl;
	s << "                     (default is " << defaultNumBits << ")" << endl;
	s << "  --uncompressed     make the filter with uncompressed bit vector(s)" << endl;
	s << "                     (this is the default)" << endl;
	s << "  --rrr              make the filter with RRR-compressed bit vector(s)" << endl;
	s << "  --roar             make the filter with roar-compressed bit vector(s)" << endl;
	s << "  --stats[=<filename>] write bloom filter stats to a text file" << endl;
	s << "                     (if no filename is given this is derived from the bloom" << endl;
	s << "                     filter filename)" << endl;
	s << endl;
	s << "When --list is used, each line of the file corresponds to a bloom filter. The" << endl;
	s << "format of each line is" << endl;
	s << "  <filename> [<filename>..] [--kmersin] [--out=<filename>]" << endl;
	s << "with meaning the same as on the command line. No other options (e.g. --k or" << endl;
	s << "--bits) are allowed in the file. These are specified on the command line and" << endl;
	s << "will affect all the bloom filters." << endl;
	s << endl;
	s << "When --kmersin is used, each line of the sequence input files is a single kmer," << endl;
	s << "as the first field in the line. Any additional fields on the line are ignored." << endl;
	s << "For example, with --k=20 this might be" << endl;
	s << "  ATGACCAGATATGTACTTGC" << endl;
	s << "  TCTGCGAACCCAGACTTGGT" << endl;
	s << "  CAAGACCTATGAGTAGAACG" << endl;
	s << "   ..." << endl;
	s << "Every kmer in the file(s) is added to the filter. No counting is performed," << endl;
	s << "and --min is not allowed." << endl;
	s << endl;
	s << "When --bits=<N> and --modulus=<M> are both used, the M bits that would normally" << endl;
	s << "be in the bloom filter are subsampled, and only the first N bits are populated." << endl;
	s << "N must be no larger than M. N can be set as a fraction or percentage of M, e.g." << endl;
	s << "--bits=1/20 or --bits=5%." << endl;
	s << endl;
	//    123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	s << "If more than one --bits=<N> and more than one --modulus=<M> are used, there" << endl;
	s << "must be the same number of each. If more than one --modulus=<M> or --bits=<N>" << endl;
	s << "is used, the number of bits and modulus will be appended to the default output" << endl;
	s << "filenames. If --out is used instead, it must contain at least one of the" << endl;
	s << "substrings {bits} or {modulus}" << endl;
	}

void MakeBFCommand::debug_help
   (std::ostream& s)
	{
	s << "--debug= options" << endl;
	s << "  settings" << endl;
	s << "  add" << endl;
	s << "  contains" << endl;
	s << "  kmers" << endl;
	s << "  strings" << endl;
	s << "  v1file" << endl;
	}

void MakeBFCommand::parse
   (int		_argc,
	char**	_argv)
	{
	int		argc;
	char**	argv;
	vector<float> numBitsFraction;

	// defaults

	listFilename  = "";
	inputIsKmers  = false;
	bfFilename    = "";
	kmerSize      = defaultKmerSize;      bool kmerSizeSet       = false;
	minAbundance  = defaultMinAbundance;       minAbundanceSet   = false;
	numThreads    = defaultNumThreads;
	numHashes     = defaultNumHashes;     bool numHashesSet      = false;
	hashSeed1     = 0;                    bool hashSeed1Set      = false;
	compressor    = bvcomp_uncompressed;  bool compressorSet     = false;
	hashSeed2     = 0;                    bool hashSeed2Set      = false;

	numFilters = 0;
	numBits.clear();
	numBitsFraction.clear();
	hashModulus.clear();

	outputStats   = false;
	statsFilename = "";

#ifdef useJellyHash
	hashSeed1 = JellyHashSeed;
#endif

	// skip command name

	argv = _argv+1;  argc = _argc - 1;
	if (argc <= 0) chastise ();

	//////////
	// scan arguments
	//////////

	for (int argIx=0 ; argIx<argc ; argIx++)
		{
		string arg = argv[argIx];
		string argVal;
		if (arg.empty()) continue;

		string::size_type argValIx = arg.find('=');
		if (argValIx == string::npos) argVal = "";
		                         else argVal = arg.substr(argValIx+1);

		// --help, etc.

		if ((arg == "--help")
		 || (arg == "-help")
		 || (arg == "--h")
		 || (arg == "-h")
		 || (arg == "?")
		 || (arg == "-?")
		 || (arg == "--?"))
			{ usage (cerr);  std::exit (EXIT_SUCCESS); }

		if ((arg == "--help=debug")
		 || (arg == "--help:debug")
		 || (arg == "?debug"))
			{ debug_help(cerr);  std::exit (EXIT_SUCCESS); }

		// --kmersin

		if ((arg == "--kmersin")
		 ||	(arg == "--askmers="))
			{ inputIsKmers = true;  continue; }

		// --out=<filename>

		if ((is_prefix_of (arg, "--out="))
		 ||	(is_prefix_of (arg, "--output=")))
			{ bfFilename = argVal;  continue; }

		// --list=<filename>

		if (is_prefix_of (arg, "--list="))
			{ listFilename = argVal;  continue; }

		// --asper=<filename>

		if (is_prefix_of (arg, "--asper="))
			{ asPerFilename = argVal;  continue; }

		// --k=<N>

		if ((is_prefix_of (arg, "K="))
		 ||	(is_prefix_of (arg, "--K="))
		 ||	(is_prefix_of (arg, "k="))
		 ||	(is_prefix_of (arg, "--k="))
		 ||	(is_prefix_of (arg, "--kmer="))
		 ||	(is_prefix_of (arg, "--kmersize=")))
			{
			kmerSize = string_to_u32(argVal);
			kmerSizeSet = true;
			continue;
			}

		// --min=<N>

		if ((is_prefix_of (arg, "--min="))
		 ||	(is_prefix_of (arg, "--abundance=")))
			{
			minAbundance = string_to_u32(argVal);
			if (minAbundance < 1) minAbundance = 1;
			minAbundanceSet = (minAbundance > 1);
			continue;
			}

		// --threads=<N>

		if ((is_prefix_of (arg, "--threads="))
		 ||	(is_prefix_of (arg, "T="))
		 ||	(is_prefix_of (arg, "--T=")))
			{
			numThreads = string_to_u32(argVal);
			if (numThreads == 0)
				chastise ("(in \"" + arg + "\") number of threads cannot be zero");
			continue;
			}

		// --hashes=<N>

		if ((is_prefix_of (arg, "--hashes="))
		 ||	(is_prefix_of (arg, "H="))
		 ||	(is_prefix_of (arg, "--H=")))
			{
			numHashes = string_to_u32(argVal);
			numHashesSet = true;
			continue;
			}

		// --seed=<number> or --seed=<number>,<number>

		if ((is_prefix_of (arg, "--seed="))
		 || (is_prefix_of (arg, "S="))
		 || (is_prefix_of (arg, "--S=")))
			{
			string::size_type commaIx = argVal.find(',');
			if (commaIx == string::npos)
				{
			    hashSeed1 = string_to_u64 (argVal);
				hashSeed1Set = true;
				}
			else
				{
			    hashSeed1 = string_to_u64 (argVal.substr(0,commaIx));
			    hashSeed2 = string_to_u64 (argVal.substr(commaIx+1));
				hashSeed1Set = hashSeed2Set = true;
				}
			continue;
			}

		// --modulus=<M>

		if ((is_prefix_of (arg, "--modulus="))
		 ||	(is_prefix_of (arg, "M="))
		 ||	(is_prefix_of (arg, "--M=")))
			{
			hashModulus.emplace_back(string_to_unitized_u64(argVal));
			continue;
			}

		// --bits=<fraction> or --bits=<percentage> or --bits=<N>

		if ((is_prefix_of (arg, "--bits="))
		 ||	(is_prefix_of (arg, "B="))
		 ||	(is_prefix_of (arg, "--B=")))
			{
			if ((is_suffix_of (arg, "%"))
			 || (arg.find_last_of("/") != string::npos))  // arg contains "/"
				{
				numBits.emplace_back(0);
				numBitsFraction.emplace_back(string_to_probability(argVal));
				}
			else
				{
				numBits.emplace_back(string_to_unitized_u64(argVal));
				numBitsFraction.emplace_back(0.0);
				}
			continue;
			}

		// bit vector type

		if (arg == "--uncompressed")
			{
			compressor    = bvcomp_uncompressed;
			compressorSet = true;
			continue;
			}

		if ((arg == "--rrr")
		 || (arg == "--RRR"))
			{
			compressor    = bvcomp_rrr;
			compressorSet = true;
			continue;
			}

		if ((arg == "--roar")
		 || (arg == "--roaring"))
			{
			compressor    = bvcomp_roar;
			compressorSet = true;
			continue;
			}

		// (unadvertised) special bit vector types

		if ((arg == "--zeros")
		 || (arg == "--allzeros")
		 || (arg == "--all_zeros")
		 || (arg == "--all-zeros"))
			{
			compressor    = bvcomp_zeros;
			compressorSet = true;
			continue;
			}

		if ((arg == "--ones")
		 || (arg == "--allones")
		 || (arg == "--all_ones")
		 || (arg == "--all-ones"))
			{
			compressor    = bvcomp_ones;
			compressorSet = true;
			continue;
			}

		if (arg == "--uncrrr")
			{
			compressor    = bvcomp_unc_rrr;
			compressorSet = true;
			continue;
			}

		if (arg == "--uncroar")
			{
			compressor    = bvcomp_unc_roar;
			compressorSet = true;
			continue;
			}

		// --stats[=<filename>]

		if (arg == "--stats")
			{ outputStats = true; continue; }

		if ((is_prefix_of (arg, "--stats=")))
			{
			outputStats   = true;
			statsFilename = argVal;
			continue;
			}

		// (unadvertised) debug options

		if (arg == "--debug")
			{ debug.insert ("debug");  continue; }

		if (is_prefix_of (arg, "--debug="))
			{
		    for (const auto& field : parse_comma_list(argVal))
				debug.insert(to_lower(field));
			continue;
			}

		// unrecognized --option

		if (is_prefix_of (arg, "--"))
			chastise ("unrecognized option: \"" + arg + "\"");

		// <filename>

		seqFilenames.emplace_back(strip_blank_ends(arg));
		}

	// if an "as per" file was given, modify settings that weren't otherwise
	// specified

	if (not asPerFilename.empty())
		{
		BloomFilter bf(asPerFilename);
		bf.load();

		if (!kmerSizeSet)    { kmerSize    = bf.kmerSize;              kmerSizeSet    = true; }
		if (!numHashesSet)   { numHashes   = bf.numHashes;             numHashesSet   = true; }
		if (!hashSeed1Set)   { hashSeed1   = bf.hashSeed1;             hashSeed1Set   = true; }
		if (!hashSeed2Set)   { hashSeed2   = bf.hashSeed2;             hashSeed2Set   = true; }
		// øøø don't access bvs like this
		if (!compressorSet)  { compressor  = bf.bvs[0]->compressor();  compressorSet  = true; }
		if (hashModulus.size() == 0)
			hashModulus.emplace_back(bf.hashModulus);
		if (numBits.size() == 0)
			{
			numBits.emplace_back(bf.numBits);
			numBitsFraction.emplace_back(0.0);
			}

		}

	// sanity checks

	if ((compressor == bvcomp_zeros)
	 || (compressor == bvcomp_ones))
		{
		if (not listFilename.empty())
			chastise ("cannot use --list with --zeros or --ones");
		if (not seqFilenames.empty())
			chastise ("cannot use sequence files with --zeros or --ones");
		if (bfFilename.empty())
			chastise ("--zeros or --ones requires --out");
		}
	else if (listFilename.empty())
		{
		if (seqFilenames.empty())
			chastise ("at least one sequence filename is required");
		}
	else
		{
		if (not seqFilenames.empty())
			chastise ("cannot use --list with sequence filenames (e.g. " + seqFilenames[0] + ") in the command");

		if (not bfFilename.empty())
			chastise ("cannot use --list with a filter filename (" + bfFilename + ") in the command");
		if (not statsFilename.empty())
			chastise ("cannot use --list with a stats filename (" + statsFilename + ") in the command");
		}

	if (kmerSize == 0)
		chastise ("kmer size cannot be zero");

	if (numHashes == 0)
		chastise ("number of hash functions cannot be zero");

	if (numHashes == 1)
		hashSeed2 = 0;
	else if (not hashSeed2Set)
		hashSeed2 = hashSeed1 + 1;

	if ((inputIsKmers) and (minAbundanceSet))
		chastise ("cannot use --kmersin with  --min");

	// resolve hashModulus and numBits

	if (numBits.size() == 0)
		{
		numBits.emplace_back((u64)defaultNumBits); // (the u64 cast placates the compiler)
		numBitsFraction.emplace_back(0.0);
		}

	if ((hashModulus.size() > 1) && (numBits.size() > 1))
		{
		if (hashModulus.size() != numBits.size())
			chastise ("can't understand "
                    + std::to_string(hashModulus.size()) + " instances of --modulus=<N> with "
                    + std::to_string(numBits.size()) + " instances of --bits=<N>");
		}

	if ((hashModulus.size() == 1) && (numBits.size() > 1))
		{
		for (size_t ix=1 ; ix<numBits.size() ; ix++)
			hashModulus.emplace_back(hashModulus[0]);
		}

	if ((hashModulus.size() > 1) && (numBits.size() == 1))
		{
		for (size_t ix=1 ; ix<hashModulus.size() ; ix++)
			{
			numBits.emplace_back(numBits[0]);
			numBitsFraction.emplace_back(numBitsFraction[0]);
			}
		}

	bool anyNumBitsIsFraction = false;
	for (size_t ix=0 ; ix<numBits.size() ; ix++)
		{
		if (numBitsFraction[ix] > 0.0)
			{ anyNumBitsIsFraction = true; break; }
		}

	if (anyNumBitsIsFraction)
		{
		if (hashModulus.size() == 0)
			chastise ("fractional bit count requires a hash modulus (--modulus=<N>)");

		for (size_t ix=0 ; ix<std::max(hashModulus.size(),numBits.size()) ; ix++)
			{
			size_t mIx = (ix < hashModulus.size())? ix : 0;
			size_t bIx = (ix < numBits.size())? ix : 0;
			if (numBitsFraction[bIx] == 0.0) continue;
			u64 nBits = ceil(numBitsFraction[bIx] * hashModulus[mIx]);
			if (nBits < 2) nBits = 2;
			numBits[bIx] = nBits;
			//numBitsFraction[bIx] = 0.0;  // (unnecessary)
			}
		}

	for (size_t ix=0 ; ix<numBits.size() ; ix++)
		{
		if (numBits[ix] < 2)
			chastise ("number of bits must be at least 2");
		}

	if (hashModulus.size() == 0)
		{
		for (size_t ix=0 ; ix<numBits.size() ; ix++)
			hashModulus.emplace_back(numBits[ix]);
		}

	for (size_t ix=0 ; ix<numBits.size() ; ix++)
		{
		if (hashModulus[ix] < numBits[ix])
			chastise ("hash modulus (" + std::to_string(hashModulus[ix]) + ")"
			        + " cannot be less than the number of bits"
			        + " (" + std::to_string(numBits[ix]) + ")");
		}

	numFilters = numBits.size();

	if ((numFilters > 1) && (not bfFilename.empty()))
		{
		string::size_type fieldIx = bfFilename.find("{bits}");
		if (fieldIx == string::npos)
			fieldIx = bfFilename.find("{modulus}");
		if (fieldIx == string::npos)
			chastise ("because more than one instance of --modulus=<N> or --bits=<N> is used,"
			          "\nthe output filename must contain {bits} and/or {modulus}");
		}

	if ((numFilters > 1) && (not statsFilename.empty()))
		chastise ("because more than one instance of --modulus=<N> or --bits=<N> is used,"
		          "\ncannot use a stats filename (" + statsFilename + ") in the command");

	// report settings

	if (contains(debug,"settings"))
		{
		cerr << "kmerSize    = " << kmerSize    << endl;
		cerr << "numHashes   = " << numHashes   << endl;
		cerr << "hashSeed1   = " << hashSeed1   << endl;
		cerr << "hashSeed2   = " << hashSeed2   << endl;
		if (numFilters == 1)
			{
			cerr << "hashModulus = " << hashModulus[0] << endl;
			cerr << "numBits     = " << numBits[0]     << endl;
			}
		else
			{
			for (size_t ix=0 ; ix<numFilters ; ix++)
				{
				cerr << "hashModulus[" << ix << "] = " << hashModulus[ix] << endl;
				cerr << "numBits[" << ix << "]     = " << numBits[ix]     << endl;
				}
			}
		cerr << "compressor  = " << compressor  << endl;
		}

	return;
	}


int MakeBFCommand::execute()
	{
	// if we're to make a single filter, do so

	if (listFilename.empty())
		{
		if (inputIsKmers) make_bloom_filter_kmers ();
		             else make_bloom_filter_fasta ();
		}

	// otherwise, make a series of filters according to each line specified in
	// a file

	else
		{
		std::ifstream in (listFilename);
	    if (not in)
			fatal ("error: failed to open \"" + listFilename + "\"");

		string line;
		int lineNum = 0;
		while (std::getline (in, line))
			{
			lineNum++;
			seqFilenames.clear();
			inputIsKmers = false;
			bfFilename   = "";

			vector<string> tokens = tokenize(line);
			for (size_t argIx=0 ; argIx<tokens.size() ; argIx++)
				{
				string arg = tokens[argIx];
				string argVal;

				string::size_type argValIx = arg.find('=');
				if (argValIx == string::npos) argVal = "";
				                         else argVal = arg.substr(argValIx+1);

				if ((arg == "--kmersin")
				 ||	(arg == "--askmers="))
					{
					if (minAbundanceSet)
						fatal ("cannot use --kmersin, with --min on the command line"
						       " (at line " + std::to_string(lineNum)
						     + " in " + listFilename + ")");
					inputIsKmers = true;
					continue;
					}

				if ((is_prefix_of (arg, "--out="))
				 ||	(is_prefix_of (arg, "--output=")))
					{ bfFilename = argVal;  continue; }

				if (is_prefix_of (arg, "--"))
					fatal ("unrecognized field: \"" + arg + "\""
					     + " at line " + std::to_string(lineNum)
					     + " in " + listFilename);

				seqFilenames.emplace_back(strip_blank_ends(arg));
				}

			if (inputIsKmers) make_bloom_filter_kmers ();
			             else make_bloom_filter_fasta ();
			}

		in.close();
		}

	return EXIT_SUCCESS;
	}


void MakeBFCommand::make_bloom_filter_fasta()  // this also supports fastq
	{
	assert (hashModulus.size() == numFilters);
	assert (numBits.size()     == numFilters);

	vector<string> bfOutFilename;
	for (size_t ix=0 ; ix<numFilters ; ix++)
		bfOutFilename.emplace_back(build_output_filename(ix));

	// create the hash table, with jellyfish defaults

	const u64 hashSize    = 10*1000*1000;
	const u32 numReprobes = 126;
	const u32 counterLen  = 7;

	unsigned int savedKmerSize = jellyfish::mer_dna::k();
	jellyfish::mer_dna::k(kmerSize);

	mer_hash_type merHash (hashSize, kmerSize*2, counterLen, numThreads, numReprobes);

	// count the kmers
	// nota bene: MerCounter internally discards kmers containing any non-ACGT
	// $$$ ERROR_CHECK need to trap exceptions from the jellyfish stuff
	// $$$ ERROR_CHECK does jellyfish give us any indication if one of the sequence files doesn't exist?

	MerCounter counter(numThreads, merHash, seqFilenames.begin(), seqFilenames.end());
	counter.exec_join (numThreads);

	// build the bloom filters

	vector<BloomFilter*> bf;
	for (size_t ix=0 ; ix<numFilters ; ix++)
		{
		bf.emplace_back(new BloomFilter(bfOutFilename[ix], kmerSize,
		                                numHashes, hashSeed1, hashSeed2,
		                                numBits[ix], hashModulus[ix]));
		if (contains(debug,"add"))      bf[ix]->dbgAdd      = true;
		if (contains(debug,"contains")) bf[ix]->dbgContains = true;

		bf[ix]->new_bits (compressor);
		}

	const auto jfAry = merHash.ary();
	const auto end   = jfAry->end();
	u64 kmersAdded = 0;
	for (auto kmer=jfAry->begin() ; kmer!=end ; ++kmer)
		{
		auto& keyValuePair = *kmer;
		if (keyValuePair.second >= minAbundance)
			{
			if (contains(debug,"kmers"))
				cerr << keyValuePair.first << " " << keyValuePair.second << endl;

			for (size_t ix=0 ; ix<numFilters ; ix++)
				{
				if (contains(debug,"strings"))
					bf[ix]->add (keyValuePair.first.to_str());
				else
					bf[ix]->add ((u64*) keyValuePair.first.data());
				}

			kmersAdded++;
			}
		}

	jellyfish::mer_dna::k(savedKmerSize);	// restore jellyfish kmer size

	for (size_t ix=0 ; ix<numFilters ; ix++)
		{
		if (not contains(debug,"v1file"))
			{
			bf[ix]->setSizeKnown = true;
			bf[ix]->setSize      = kmersAdded;
			}

		if ((compressor == bvcomp_unc_rrr)
		 || (compressor == bvcomp_unc_roar))
			{
			BitVector* bv = bf[ix]->bvs[0];
			bv->unfinished();
			}

		bf[ix]->reportSave = true;
		bf[ix]->save();
		delete bf[ix];
		}

	// report stats to a file and/or the console

	if ((outputStats) or (contains(debug,"fprate")))
		{
		for (size_t ix=0 ; ix<numFilters ; ix++)
			report_stats(bfOutFilename[ix],kmersAdded,ix);
		}
	}


void MakeBFCommand::make_bloom_filter_kmers()
	{
	assert (hashModulus.size() == numFilters);
	assert (numBits.size()     == numFilters);

	vector<string> bfOutFilename;
	for (size_t ix=0 ; ix<numFilters ; ix++)
		bfOutFilename.emplace_back(build_output_filename(ix));

	// build the bloom filters

	vector<BloomFilter*> bf;
	for (size_t ix=0 ; ix<numFilters ; ix++)
		{
		bf.emplace_back(new BloomFilter(bfOutFilename[ix], kmerSize,
		                                numHashes, hashSeed1, hashSeed2,
		                                numBits[ix], hashModulus[ix]));
		if (contains(debug,"add"))      bf[ix]->dbgAdd      = true;
		if (contains(debug,"contains")) bf[ix]->dbgContains = true;

		bf[ix]->new_bits (compressor);
		}

	u64 kmersAdded = 0;
	for (const auto& kmersFilename : seqFilenames)
		{
		std::ifstream in (kmersFilename);
		if (not in)
			fatal ("error: failed to open \"" + kmersFilename + "\"");

		string line, kmer;
		int lineNum = 0;
		while (std::getline (in, line))
			{
			lineNum++;
			line = strip_blank_ends (line);
			if (line.empty()) continue;

			string::size_type whitespaceIx = line.find_first_of(" \t");
			if (whitespaceIx == string::npos) kmer = line;
			                             else kmer = line.substr(0,whitespaceIx);

			if (kmer.length() != kmerSize)
				fatal ("error: expected " + std::to_string(kmerSize) + "-mer"
				     + " but encountered " + std::to_string(kmer.length()) + "-mer"
				     + " (at line " + std::to_string(lineNum)
				     + " in " + kmersFilename + ")");

			string::size_type badIx = kmer.find_first_not_of("ACGTacgt");
			if (badIx != string::npos) continue;

			if (contains(debug,"kmers"))
				cerr << kmer << endl;

			for (size_t ix=0 ; ix<numFilters ; ix++)
				bf[ix]->add (kmer);
			kmersAdded++;
			}
		}

	for (size_t ix=0 ; ix<numFilters ; ix++)
		{
		if (not contains(debug,"v1file"))
			{
			bf[ix]->setSizeKnown = true;
			bf[ix]->setSize      = kmersAdded;
			}

		if ((compressor == bvcomp_unc_rrr)
		 || (compressor == bvcomp_unc_roar))
			{
			BitVector* bv = bf[ix]->bvs[0];
			bv->unfinished();
			}

		bf[ix]->reportSave = true;
		bf[ix]->save();
		delete bf[ix];
		}

	// report stats to a file and/or the console

	if ((outputStats) or (contains(debug,"fprate")))
		{
		for (size_t ix=0 ; ix<numFilters ; ix++)
			report_stats(bfOutFilename[ix],kmersAdded,ix);
		}
	}


void MakeBFCommand::report_stats(const string& bfOutFilename, u64 kmersAdded,size_t filterNum)
	{
	double fpRate = BloomFilter::false_positive_rate(numHashes,hashModulus[filterNum],kmersAdded);

	if (outputStats)
		{
		string statsOutFilename = build_stats_filename(bfOutFilename);
		cerr << "writing bloom filter stats to \"" << statsOutFilename << "\"" << endl;

	    std::ofstream statsF(statsOutFilename);

		statsF << "#filename"
		       << "\tnumHashes"
		       << "\tnumBits"
		       << "\tkmersAdded"
		       << "\tbfFpRate"
		       << endl;
		statsF <<         bfOutFilename
		       << "\t" << numHashes
		       << "\t" << hashModulus[filterNum]
		       << "\t" << kmersAdded
		       << "\t" << fpRate
		       << endl;
		}

	if (contains(debug,"fprate"))
		{
		cerr << bfOutFilename << " kmers inserted: " << kmersAdded << endl;
		cerr << bfOutFilename << " estimated BF false positive rate: " << fpRate << endl;
		}
	}


string MakeBFCommand::build_output_filename(size_t filterNum)
	{
	string bfOutFilename = bfFilename;

	if (bfOutFilename.empty())
		{
		string ext = "." + BitVector::compressor_to_string(compressor) + ".bf";
		if (ext == ".uncompressed.bf") ext = ".bf";

		if (numFilters > 1)
			{
			string bfSizeText = std::to_string(numBits[filterNum]) + "_of_" + std::to_string(hashModulus[filterNum]);
			ext = "." + bfSizeText + ext;
			}

		string seqFilename = seqFilenames[0];
		string::size_type dotIx = seqFilename.find_last_of(".");
		if (dotIx == string::npos)
			bfOutFilename = seqFilename + ext;
		else
			bfOutFilename = seqFilename.substr(0,dotIx) + ext;
		}
	else if (numFilters > 1)
		{
		// note that command line parsing required that bfFilename contain
		// at least one of {bits} or {modulus} when numFilters>1
		string::size_type fieldIx = bfOutFilename.find("{bits}");
		if (fieldIx != string::npos)
			bfOutFilename = bfOutFilename.replace(fieldIx,6,std::to_string(numBits[filterNum]));
		fieldIx = bfOutFilename.find("{modulus}");
		if (fieldIx != string::npos)
			bfOutFilename = bfOutFilename.replace(fieldIx,9,std::to_string(hashModulus[filterNum]));
		}

	return bfOutFilename;
	}


string MakeBFCommand::build_stats_filename(const string& bfOutFilename)
	{
	string statsOutFilename = statsFilename;

	if (statsOutFilename.empty())
		{
		string ext = ".stats";

		string::size_type dotIx = bfOutFilename.find_last_of(".");
		if (dotIx == string::npos)
			statsOutFilename = bfOutFilename + ext;
		else
			statsOutFilename = bfOutFilename.substr(0,dotIx) + ext;
		}

	return statsOutFilename;
	}

