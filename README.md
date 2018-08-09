Strain-GeMS
===========

Strain-GeMS: identifying strains within metagenomic species

One line pitcher
===========
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

If you already have MetaPhlAn results for some(or, all) of the samples, you can specify them to save some time as:

    //
    sample: [sample1_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    metaphlan2: [metaphlan2 results for sample1]
    //
    sample: [sample2_ID]
    fq1: [forward reads fastq]
    fq2: [reverse/mate reads fastq]
    metaphlan2: [metaphlan2 results for sample2]


 
Below is a detailed usage of Strain-GeMS.py:
  
  Options:
  
    --version             show program's version number and exit
  
    -h, --help            show this help message and exit

  Compulsory parameters:
    
    There options are compulsory, and may be supplied in any order.

    -o DIR, --outdir=DIR
                        The output directory of Strain-GeMS.
    -c FILE, --config=FILE
                        The configuration file of the Strain-GeMS project.

  Optional parameters:
  
    There options are optional, and may be supplied in any order.

    -t INT, --num_proc=INT
                        Number of processor for Strain-GeMS to use [default: 1].
    -d STRING, --ref_db=STRING
                        The prefix of species reference. [default:
                        Strain-GeMS/db/ref_db].
    -g STRING, --gsize_db=STRING
                        The directory of species average genome size DB.
                        [default: Strain-GeMS/db/gsize.db].
    --bowtie2=STRING    Path to bowtie2 binary, specify if not in env path
                        [default: bowtie2].
                        Bowtie2 citation: Langmead B. and Salzberg S., Nat.
                        Methods, 2012.
                        Bowtie2 page: http://bowtie-
                        bio.sourceforge.net/bowtie2
    --bowtie2_build=STRING
                        Path to bowtie2-build bin  ary, specify if not in env
                        path [default: bowtie2-build].
                        Bowtie2 citation: Langmead B. and Salzberg S., Nat.
                        Methods, 2012.
                        Bowtie2 page: http://bowtie-
                        bio.sourceforge.net/bowtie2
    --samtools=STRING   Path to samtools binary, specify if not in env path
                        [default: samtools].
                        Samtools citation: Li H., et al, Bioinformatics, 2009.
                        Samtools webpage: http://samtools.sourceforge.net
    -m STRING, --metaphlan2=STRING
                        Path to metaphlan2 script, specify if not in env path
                        [default: metaphlan.py].
                        MetaPhlAn citation: Segata N. et al, Nat. Methods,
                        2012.
                        MetaPhlAn2 page:
                        http://huttenhower.sph.harvard.edu/metaphlan2

  Straining parameters:
  
    There options are optional, and may be supplied in any order.
    
    --min_cov=FLOAT     Minimum coverage of a species in a sample to be
                        considered [default: 10, range: 5+].

  Runtime settings:
  
    There options are optional, and may be supplied in any order.

    -q, --quiet         Suppress printing detailed runtime information,
                        only important warnings/errors will show [default:
                        False].



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
 


 
