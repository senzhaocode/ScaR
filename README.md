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
        
# 2. Usage of the "select_read.pl"
  2.1 Look at running parameters
      
      perl select_read.pl --help

  2.2 An example of running:
      
      perl select_read.pl \
      
      --first ~/examples/input/raw_1.fastq \ 
      # fastq file path for the first end of paired-end reads
      
      --second ~/examples/input/raw_2.fastq \
      # fastq file path for the second end of paired-end reads
      
      --geneA RCC1 --geneB ABHD12B \
      # fusion partner gene names (Refseq gene symbol and Ensembl id are accepted, but never mix them together)
      
      --anchor 6 \ 
      # (default: 6)
      # Set the length of anchor. The minimum number of bases is required to match to geneA/geneB region in the scaffold sequence.
      # If users want to more specificity of read mapping, just increase this value.
      
      --trimm 0 \
      # (default: 0)
      # Set whether raw input fastq reads are trimmed (1) or not (0)
      
      --length 48 \
      # Set the maximum length of fastq read, this option is only available when raw fastq reads are trimmed (--trimm 1)
      
      --p 8 \ 
      # (default: 8)
      # The number of threads, and make sure that it should be the same as the number of CPUs allocated in jobscript
      
      --input ~/data/ \
      # Set input path of genomic data and annotation
      
      --output ~/examples/output/ \
      # Output directory
      
      --scaffold ~/examples/input/RCC1_ABHD12B_scaff_seq.fa \
      # Fusion scaffold sequences, users can extract them from raw output files of de novo fusion finders (e.g. deFuse, fusioncatcher, SOAPfuse). A list of candidate sequences in fasta format are accepted. 
      # For instance:
      #	>alt_0
      #	XXXXXXXXXXXXXXXXXX|YYYYYYYYYYYYYYYY
      #	>alt_1
      #	XXXXXXXXXXXXXXXXXX*YYYYYYYYYYYYYYY
      # NOTE: 1. It only accepts '*' or '|' as a separator for breakpoint sequences. Sequences 'XXXXXXXXXXX' and 'YYYYYYYYY' should correspond to geneA and geneB, respectively.
      #       2. In order to ensure the specificity of breakpoint sequence, we recommend that 'XXXXXXXX' and 'YYYYYYYY' should be at least 20 bp.
      #       3. In current version, the breakpoint sequence should be composed of cDNA (i.e. exon region). If it contains intron/intergenic sequences, the program will stop running scaffold realignment.

  2.3 Some tips for run scripts using slurm jobscript.
  
      Generally, "select_read.pl" script does not need huge memory, instead the implementation of multiple number of threads will speed up running process. For example: the setting "SBATCH --cpus-per-task=8 and SBATCH --mem-per-cpu=1G" should be much more efficient than the setting "SBATCH --cpus-per-task=4 and SBATCH --mem-per-cpu=4G"
      
# 3. Output results
  For example: ~/examples/output/
    
    |--- "scaffold_*_seq.fa" (cDNA sequences of geneA and geneB, and breakpoint sequence of scaffold in fasta format)
    |--- "hisats_noclip.sorted.bam" (BAM file: done by Hisat2 no-splicing alignment model)
    |--- "discordant_split_1.txt, discordant_split_2.txt" (paired-end reads in fastq format: one-end maps to scaffold; the other maps to cDNA sequences of geneA/geneB -- extracted from hisats_noclip.sorted.bam)
    |--- "singlton_split_1.txt, singlton_split_2.txt" (paired-end reads in fastq format: one-end maps to scaffold; the other shows no mapping to cDNA sequences of geneA/geneB -- extracted from hisats_noclip.sorted.bam)
    |--- "spanning_1.txt, spanning_2.txt" (paired-end reads in fastq format: one-end maps to cDNA sequence of geneA; the other maps to cDNA sequence of geneB -- extracted from hisats_noclip.sorted.bam)
    |--- "read_mapped_info" (mapping summary of discordant/singlton split reads and spanning reads)
    |
   *|--- "final_read_mapped_info" (mapping summary of filtered discordant/singlton split reads and filtered spanning reads)
   *|--- "final_split_1.txt, final_split_2.txt" (paired-end reads in fastq format: combine discordant and singlton split reads after filtering out unspecific read mapping at the genome level)
    |--- "final_spanning_1.txt, final_spanning_1.txt" (paired-end reads in fastq format: spanning reads after filtering out unspecific mapping at the genome level)
   *|--- "final_noclip.sorted.bam" (show filtered discordant and singlton split reads, and spanning reads mapped to scaffold sequence, done by hisat2 no-splicing alignment model)
    |    If there are no spanning and discordant/singlton split reads, the "final_noclip.sorted.bam" is not present.
    |--- "tmp" folder (if users want to look at more detail of processing steps):
        |--- "hisats_noclip.sam" (SAM format of "hisats_noclip.sorted.bam")
        |--- "spanning_sec.sam" (align "spanning_1.txt, spanning_2.txt" to genome reference, done by hisat2 splicing alignment model)
        |--- "discordant_split_sec.sam" (align "discordant_split_1.txt, discordant_split_2.txt" to genome reference, done by hisat2 splicing alignment model)
        |--- "singlton_split_sec.sam" (align "singlton_split_1.txt, singlton_split_2.txt" to genome reference, done by hisat2 splicing alignment model)
    *  most important output files for users

