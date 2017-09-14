# scaffold_realign package
Use scaffold realigning strategy to detect the recurrence of a list fusion transcripts across samples

# 1. Requirements
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
    
  1.6 Genomic data and annotations
  
      ~/data/Gene_hg38.txt: gene annotation file (This file has been present in ~/data when you downlod the package)
      ~/data/transcript_cdna.fa: user have to download the "transcript_cdna.fa" (transcriptome sequence) in ~/data following the command
        cd ~/data
        wget "http://folk.uio.no/senz/Transcript_cdna.fa"
        
