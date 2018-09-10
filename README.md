Strain-GeMS
===========

Strain-GeMS: identifying strains within metagenomic species

Strain-GeMS aims to profile the strains within a species using short read metagenomes.
It deploys SNP patterns based on statistically sound SNP calling and optimzied SNP clustering 
to recover strain genotypes within a species using a set of genes.




Install
===========

You can download the zip archive of it from : 

https:https://github.com/HUST-NingKang-Lab/straingems
	
Dependencies
============

* Python-2.7 or above

+ Python libraries:

>BioPython

>Numpy 

>Scipy

>NetworkX 1.11
  
+ Third party pipelines:  

>Bowtie2 2.2.1+

>Metaphlan2 2.1+

>multigems 0.1+

>samtools 0.1.19

Note: older versions of these dependencies might work but it's not guaranteed.



Usage
===========

The basic Strain-GeMS analysis runs as below:

    Strain-GeMS.py [options] -c/--conf <config.file> -o/--outdir <output directory>

The format config file follows:

    //
    sample: [sample1_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    //
    sample: [sample2_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    ...
 
Or, if only one fastq file per sample, the config file should look like this:
    
    //
    sample: [sample1_ID]
    fq: [sample reads fastq]
    //
    sample: [sample2_ID]
    fq: [sample reads fastq]
    ...

 
Below is a detailed usage of Strain-GeMS.py:
  
  Options:
  
    -h, --help            show this help message and exit

  Compulsory parameters:
    
    There options are compulsory, and may be supplied in any order.

    -o DIR, --outdir=DIR
                        The output directory of Strain-GeMS.
    -c FILE, --config=FILE
                        The configuration file of the Strain-GeMS project.

  Optional parameters:
  
    -t INT, --num_proc=INT
                        Number of processor for Strain-GeMS to use [default: 1].



Interpret output
===========

In the project/output directory, you will find a folder called "results", in which you can find the below files and directories:

    Intra_sp_rel_ab.profiles     # this is a tabular file with strain relative abundance within species. See header for details.
    Overall_rel_ab.profiles      # this is a tabular file with strain relative abundance in overall samples. See header for details.
    uniGcode                     # this is the directory for all ".uniGcode" files, which are genotypes for strains.
    
The formats for "Intra_sp_rel_ab.profiles" and "Overall_rel_ab.profiles" are similar. An example is as below:
    
    # Species	strain_ID	masked_samples	sample_1   sample_2
    Escherichia_coli	str-1	NA	53.252835   37.212245
    Escherichia_coli	str-2	NA	46.747165   62.787755
    Salmonella_typhi    str-1   1   15.492194   41.212442
    Salmonella_typhi    str-2   1   38.313441   21.483291
    Salmonella_typhi    str-3   1   46.194365   37.304267
    
This means there are two species that passed the minimum coverage requirement for strain inference, and the relative abundance of each strain (E.coli has 2 strains and S.typhi has 3) are listed in sample_1 and sample_2 columns.

Strain-GeMS tries to infer strain compositions in samples with lower coverage, but these results might not be reliable. We offer a column called "masked_samples" to label such samples; the index (1-offset) of samples are comma-delimited.
 


 
 
Generation of test data
===========

Test data are generated using GemSIM_v1.6 by the GemReads.py command:

      -R Only for metagenome projects. Directory containing references.
      -a Only for metagenome projects. Species-abundance file.
      -n number of reads to produce. For paired end reads, number of pairs.
      -l length of reads. Integer value, or -l d for empirical distribution.
      -m error model file *_single.gzip or *_paired.gzip.
      -c use this flag if you wish to draw reads from a circular genome.
      -q quality score offset. Usually 33 or 64 (see manual).
      -o output file name prefix.
      -u Mean fragment length for paired end reads. -u d for empirical.
      -s standard deviation for fragment length. Use only with -u and -p.
      -p use only to create paired end reads.
	  
e.g. GemReads.py  -R ./ref/ -a abundance.txt  -n 5000000 -l d -p -u d -m ~/tools/GemSIM_v1.6/models/ill100v5_p.gzip -c -q 33 -o meta

Coverage information of the test data are given by MetaPhlAn2.


Link to test data
===========

The test data are upload to the NCBI website:
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA480741
 
