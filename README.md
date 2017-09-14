# scaffold_realign package
Use scaffold realigning strategy to detect the recurrence of a list fusion transcripts across samples

# 1. Requirements
  1.1 Perl version >= 5.10.2
  1.2 HISAT2 version 2.1.0 (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip)
      The executable files have been integrated in ~/bin/hisat2-2.1.0/, users do not have to download and install it independently.
      Users only need to add their paths to linux environment variables before running: 
        PATH=$PATH:/where_is_path/scaffold_map/bin/hisat2-2.0.5/
        export PATH
