# scaffold_realign package
Use scaffold realigning strategy to detect the recurrence of a list fusion transcripts across samples

# 1. Requirements (before running program)
  1.1 Perl version >= 5.10.2
  
  1.2 HISAT2 version 2.1.0 (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip)
  
      The executable files have been integrated in ~/bin/hisat2-2.1.0/, users do not have to download and install it independently.
      Users only need to add its path to linux environment variables before running: 
        PATH=$PATH:/where_is_path/scaffold_map/bin/hisat2-2.0.5/
        export PATH

  1.3 HISAT2 aligner HGFM index (genome reference plus transcripts based on Ensembl GRCh38 version) 
      
      Users have to download them in ~/data following the command:
        cd ~/data
        wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz"
        tar -vxf grch38_tran.tar.gz
        mv grch38_tran/genome_tran.* .
  
  1.4 Samtools version >= 1.1 (https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2)
      
      If samtools has not been installed in the server system, users have to download, compile and install it locally.
      Then, users add its path to linux environment variables before running:
        PATH=$PATH:/where_is_path_samtools
        export PATH
 
  1.5 R version >= 3.0.3 (https://cran.r-project.org)
      
      If R has not been installed in the server system, users have to donwload and install it locally.
      Then, users add it path to linux environment variables before running:
        PATH=$PATH:/where_is_path/R
        export PATH
    
  1.6 Genomic data and annotations
  
      ~/data/Gene_hg38.txt: gene annotation file (This file has been present in ~/data when you downlod the package)
      ~/data/transcript_cdna.fa: user have to download the "transcript_cdna.fa" (transcriptome sequence) in ~/data following the command
        cd ~/data
        wget "http://folk.uio.no/senz/Transcript_cdna.fa"
        
  1.7 Set the path of Perl library to environment variables
  
      PERL5LIB="$PERL5LIB:/where_is_path/scaffold_map/lib"
      export PERL5LIB
      NOTE: we recommend that users add the path of perl library to .bashrc, then "source .bashrc"
        
# 2. Run the "select_read.pl"
  2.1 Look at running parameters
      
      perl select_read.pl --help

  2.2 An example of running:
      
      perl select_read.pl \
      --first /work/projects/nn9313k/Andreas/RNA-seq/TCGA_testicular_53_102/4e42785e-6633-407d-afc6-848710e6f34e/*_1.fastq \ 
      # fastq file path for the first end of paired-end reads
      
      --second /work/projects/nn9313k/Andreas/RNA-seq/TCGA_testicular_53_102/4e42785e-6633-407d-afc6-848710e6f34e/*_2.fastq \
      # fastq file path for the second end of paired-end reads
      
      --geneA TMPRSS2 --geneB ERG \
      # fusion partner gene names (Refseq gene symbol and Ensembl id are accepted, but never mix them together)
      
      --anchor 6 \ (default: 6)
      # Set the length of anchor. The minimum number of bases is required to match to geneA/geneB region in the scaffold sequence.
      # If users want to more specificity of read mapping, just increase this value.
      
      --p 8 \ (default: 8)
      # The number of threads, and make sure that it should be the same as the number of CPUs allocated in jobscript
      
      
      --scaffold /work/projects/nn9313k/TCGA_prad/scaffold_map/TMPRSS2_ERG_scaff_seq.fa \
      # Fusion scaffold sequences, users can extract them from raw output files of de novo fusion finders (e.g. deFuse, fusioncatcher, SOAPfuse). A list of candidate sequences in fasta format are accepted. 
      # For instance:
      #	>alt_0
      #	XXXXXXXXXXXXXXXXXX|YYYYYYYYYYYYYYYY
      #	>alt_1
      #	XXXXXXXXXXXXXXXXXX*YYYYYYYYYYYYYYY
      # NOTE: 
