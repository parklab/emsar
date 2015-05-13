### README for EMSAR v1.2j
 Inquiries can be written to Soo Lee (duplexa@gmail.com).
 
 
### Table of contents 

1. Installation
2. Usage
 1. Building an rsh index
 2. Running on RNA-seq data
 3. Output files

***************************


# 1. Installation
 Extract the .tar.gz file. Then, EMSAR.v1.2j directory will be created. 
 
    cd ENSAR.v1.2j
    make

 Then, if you want, include your EMSAR.v1.2j directory to your path.


# 2. Usage
## 1) Building an rsh index

  An rsh index is a pre-built transcriptome index that EMSAR uses for fast computation. It contains the information about sequence sharing among transcripts. An rsh index is determined by transcriptome fasta file, read length (or read length range in case of SE), fragment length range (in case of PE), library type (eg. PE/SE and unstranded/stranded). If you set read length range or fragment length range to be large, it takes longer to compute the rsh index, but may be more reusable for a wide range of RNA-seq samples.

    Usage : emsar-build <options> fastafile readlength(range) outdir outprefix
    ex (SE, unstranded, readlength 76 bp) : emsar-build Hsa.GRCh37.65.gtf.fa 76 rsh human.rna.fna.SE.l76
    ex (SE, forward-stranded, readlength 50-60 bp) : emsar-build -s ssf Hsa.GRCh37.65.gtf.fa 50-60 rsh human.rna.fna.SE.l50-60.ssf
    ex (PE, unstranded, readlength 101 bp, fraglength range 1-400 bp, use 4 threads) : emsar-build --PE -p 8 Hsa.GRCh37.65.gtf.fa 101 rsh human.rna.fna.PE.l101.F1-400
    ex (PE, unstranded, readlength 101 bp, fraglength range 1-500 bp, Refseq format, use 8 threads) : emsar-build --PE -h R -p 8 -F 500 human.rna.fna 101 rsh human.rna.fna.PE.l101.F1-500

### Input files

	   fasta file : A transcriptome sequence library file in fasta format. The header can either be in Refseq format('>xx|xx|xx|name|xx') in which name or xx does not contain '|' letter, or in Ensembl format '>name' or '>name xx' in which name does not contain whitespace(space or tab) and xx is separated by whitespace from the name. Name must be a unique identifier of a sequence.

### Options
	  -P, --PE : paired-end data (default : single-end)
	  -s, --strand_type <strand_type> : set strand type ('ns','ssf','ssr' for single-end, 'ns','ssfr','ssrf' for paired-end). (default: ns(unstranded))
	  -p, --maxthread <num_thread>: number of threads to use simultaneously. Using multiple threads comes with a slight increase in memory usage. We recommend -p 4 for an optimal performance. (default : 1)
	  -F, --maxfraglen <Max_fraglen> : Maximum fragment length. Use a number that is safely large, unless you want to apply fragment length filtering. Default 400. Not applicable for SE.
	  -f, --minfraglen <Min_fraglen> : Minimum fragment length. Use a number that is safely small, unless you want to apply fragment length filtering. Default 1. Not applicable for SE
	  -h, --header <E|R> : fasta header option. E : Ensembl header(default), R : RefSeq header.
	  -k, --max_repeat <max_repeat> : the maximum number of alignments per read allowed. Reads exceeding this number are discarded.
	  -T, --print_sfa : print out (mate1) suffix array.
	  -v --verbose : make a (very) verbose output log.
	  -q --no_verbose : turn-off verbosity.

### Advanced Options
	  -b, --binsize <binsize> : binsize for indexing transcriptome (default : 5000). This bin size can be made smaller to improve speed at the cost of more memory, or vise versa.
	  -t, --taglen <taglen> : length of short sequence tags to use for constructing suffix array on a subset of substrings only. This affects only speed and memory usage. Currently, three values are supported (1,2,3). (default 2)

### Output files
    (1) prefix.rsh : an rsh index file
    (2) (with -T option) prefix.sfa : a suffix array file (mostly for developers)
    
    (prefix = "outdir/outprefix" as specified in the usage).




## 2) Running on RNA-seq data

  The actual transcript quantification can be done either with the pre-built rsh file (-I option) or a transcriptome fasta file (-x option). Note that the former is much faster. An rsh index file can also be built as a byproduct of this run by using an -R option. In this case, read length (range) will be learned directly from the alignment files (eg. bam / sam / bowtie1 output). Either a single alignment file or a list of alignment files (with -M option) can be provided as the last argument.

    Usage : emsar <options> -x fastafile outdir outprefix alignmentfile|alignmentfilelist
    Usage2 : emsar <options> -I rshfile outdir outprefix alignmentfile|alignmentfilelist
    Usage3 : bowtie command | emsar <options> [-x fastafile][-I rshfile] outdir outprefix

    ex : emsar -p 4 -h R -x human.rna.fna RNAseq sample22 sample22.bowtieout
    ex2 : emsar -p 4 -B --PE -R -x human.rna.fna RNAseq sample22 sample22.BAM
    ex3: emsar -M -p 16 -B -I human.rna.rsh RNAseq samples samples.BAMlist   ## BAM list file with -M and -B options
    ex4 : bowtie -v 2 -a -m 100 -p 4 human.rna sample22.fastq | emsar human.rna.fna RNAseq sample22

###	Input files

	   *Either a fasta file (-x) or an rsh file (-I) must be provided.
	   fasta file : A transcriptome sequence library file in fasta format. The header can either be in Refseq format('>xx|xx|xx|name|xx') in which name or xx does not contain '|' letter, or in Ensembl format '>name' or '>name xx' in which name does not contain whitespace(space or tab) and xx is separated by whitespace from the name. Name must be a unique identifier of a sequence.
	   rsh file : This file can be constructed by using -R option combined with -x fasta file. Once it is constructed, it can be used for other RNA-seq samples with the same read length for speedy calculation.

	   *Either an alignment file or a file listing alignment files must be provided.
	   alignmentfile : either a default bowtie output format (bowtieoutfile) or SAM/BAM format. SAM and BAM files must be used with -S and -B options, respectively. The SAM/BAM file should be either produced with the default sorting by bowtie or sorted by qname. We strongly recommend that bowtie is run with options -a -v 2 -m 100, without --best or --strata, without allowing any indels, for best results with EMSAR. The bowtie/SAM/BAM files can be streamed directly to EMSAR through a pipe. If mismatches were allowed, it is highly recommended that the SAM/BAM files contain the auxiliary MD flags that contains the match/mismatch information. Bowtie by default produces SAM/BAM files with the MD field.
	   alignmentfilelist : a text file containing one alignment file name per line. The alignment files can be one of the formats specified above. When a list file is used, -M (--multisample) option must be used along with the -B, -S options. All the files must be of the same read length and library type. Paired and single-end samples cannot be combined.

###	Options
	  -M, --multisample : multisample (default : single sample). When this option is used, a file containing a list of alignment files must be specified as input instead of an alignment file. The resulting expression values may be slightly different between multi-sample and single-sample runs, because a pooled fragment length distribution is used for a multi-sample run.
	  -P, --PE : paired-end data (default : single-end)
	  -s, --strand_type <strand_type> : set strand type ('ns','ssf','ssr' for single-end, 'ns','ssfr','ssrf' for paired-end). (default: ns(unstranded))
	  -S, --SAM : input file format is SAM (by default, default bowtie output)
	  -B, --BAM : input file format is BAM (by default, default bowtie output)
	  -p, --maxthread <num_thread>: number of threads to use simultaneously. Using multiple threads comes with a slight increase in memory usage. We recommend -p 4 for an optimal performance. (default : 1)
	  -F, --maxfraglen <Max_fraglen> : Maximum fragment length. Use a number that is safely large, unless you want to apply fragment length filtering. Default 400. Not applicable for SE.
	  -f, --minfraglen <Min_fraglen> : Minimum fragment length. Use a number that is safely small, unless you want to apply fragment length filtering. Default 1. Not applicable for SE.
	  -h, --header <E|R> : fasta header option. E : Ensembl header(default), R : RefSeq header.
	  -k, --max_repeat <max_repeat> : the maximum number of alignments per read allowed. Reads exceeding this number are discarded.
	  -g, --print_segments : print out segment information.
	  -T, --print_sfa : print out (mate1) suffix array.
	  -R, --print_rsh : create an rsh file from the current fasta file (must be combined with -x not -I) and the bam file (that determines read length). This rsh file can be reused for later samples for fast calculation.
	  -v --verbose : make a (very) verbose output log.
	  -q --no_verbose : turn-off verbosity.

### Advanced Options
	  -b, --binsize <binsize> : binsize for indexing transcriptome (default : 5000). This bin size can be made smaller to improve speed at the cost of more memory, or vise versa.
	  -t, --taglen <taglen> : length of short sequence tags to use for constructing suffix array on a subset of substrings only. This affects only speed and memory usage. Currently, three values are supported (1,2,3). (default 2)
	  -n, --nround <num_rounds> : number of MLE runs for computing mean and sd of FPKM. Default 4.
	  -e, --epsilon <epsilon> : epsilon value to check convergence of MLE (default : 1E-9)
	  -r, --precision <precision> : estimate precision (epsilon for step size) (default : 1E-15)
	  -i, --max_niter_mle <max_niter_mle> : maximum number of iteration before determining not converging and reinitializing for MLE (default : 10000)
	  -d, --delta <delta> : delta offset for MLE (default : 0)

### Output files
    (1) prefix.i.fpkm : abundance estimates in FPKM along with inferred read counts and TPM(transcripts per million) which is the estimated molecular fraction times 1E6.
    (2) prefix.i.fraglength_effect : fragment length effect (for your information)
    (3) (with -g option) prefix.i.segments : contains read counts and lengths ('EUMA') for all segments (for your information)
    (4) (with -R option) prefix.rsh : an rsh index file
    (5) (with -T option) prefix.sfa : a suffix array file (mostly for developers)

    (prefix = "outdir/outprefix" as specified in the usage; i = alignment file number (0-based))
