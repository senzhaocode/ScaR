ScaR package
========================

Use scaffold realigning approach to detect the prevalence and recurrence of known fusion transcripts across samples

![alt text](http://folk.uio.no/senz/Demo.png)

## 1. Requirements (before running program)
  
  1.1 Perl version >= 5.10.2
  
  1.2 HISAT2 version 2.1.0 (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip)
  
      The executable files have been integrated in ~/bin/hisat2-2.1.0/, users do not have to download and install it independently.
      The only thing is to add its path to linux environment variables before running: 
        PATH=$PATH:/where_is_path/ScaR/bin/hisat2-2.1.0/
        export PATH

  1.3 HISAT2 aligner HGFM index (genome reference plus transcripts based on Ensembl GRCh38 version) 
      
      Users have to download them in ~/reference following the command:
        cd ~/reference
        wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz"
        tar -vxf grch38_tran.tar.gz
        mv grch38_tran/genome_tran.* .
  
  1.4 Samtools version >= 1.3 (https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2)
      
      If samtools has not been installed in the server system, users have to download, compile and install it locally.
      Then, add its path to linux environment variables before running:
        PATH=$PATH:/where_is_path_samtools
        export PATH
 
  1.5 R version >= 3.0.3 (https://cran.r-project.org)
      
      If R has not been installed in the server system, users have to donwload and install it locally.
      Then, add it path to linux environment variables before running:
        PATH=$PATH:/where_is_path/R
        export PATH
    
  1.6 Genomic data and annotations
  
      ~/reference/Gene_hg38.txt: gene annotation file
      ~/reference/ensembl_transcript.fa: the human transcriptome sequences annotated from ensembl database (Ensembl Archive Release 89)
      ~/reference/gencode_transcript.fa: the human transcriptome sequences annotated from GENCODE database (Release version 27)
      ~/reference/ucsc_transcript.fa: the human transcriptome sequences annotated from UCSC database (Release date: Jan 2018)
      Typing the command
        cd ~/reference
        wget "http://folk.uio.no/senz/Gene_hg38.txt"
        wget "http://folk.uio.no/senz/ensembl_transcript.fa"
        wget "http://folk.uio.no/senz/gencode_transcript.fa"
        wget "http://folk.uio.no/senz/ucsc_transcript.fa"
        wget "http://folk.uio.no/senz/ucsc_refGene.txt"
        
  1.7 Set the path of Perl library to environment variables
  
      PERL5LIB="$PERL5LIB:/where_is_path/ScaR/lib"
      export PERL5LIB
      NOTE: we recommend that users add the path of perl libraries to .bashrc, then "source .bashrc"
        
## 2. Usage of the "select_read.pl"
  2.1 See running parameters
      
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
      
      --trans_ref ensembl
      # (default: ensembl)
      # The setting of annotation resources, users could choose others (e.g. "gencode" or "ucsc").
      
      --p 8 \ 
      # (default: 8)
      # The number of threads, and make sure that it should be the same as the number of CPUs allocated in jobscript
      
      --input ~/reference/ \
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
      
## 3. Output results
  For example: ~/examples/output/
  
  * `scaffold_*_seq.fa` (cDNA sequences of geneA and geneB, and breakpoint sequence of scaffold in fasta format)
  * `hisats_noclip.sorted.bam` (BAM file: done by Hisat2 no-splicing alignment model)
  * `discordant_split_1.txt, discordant_split_2.txt` (paired-end reads in fastq format: one-end maps to scaffold; the other maps to cDNA sequences of geneA/geneB -- extracted from hisats_noclip.sorted.bam)
  * `singlton_split_1.txt, singlton_split_2.txt` (paired-end reads in fastq format: one-end maps to scaffold; the other shows no mapping to cDNA sequences of geneA/geneB -- extracted from hisats_noclip.sorted.bam)
  * `spanning_1.txt, spanning_2.txt` (paired-end reads in fastq format: one-end maps to cDNA sequence of geneA; the other maps to cDNA sequence of geneB -- extracted from hisats_noclip.sorted.bam)
  * `read_mapped_info` (mapping summary of discordant/singlton split reads and spanning reads)
  * **`*final_read_mapped_info`** (mapping summary of filtered discordant/singlton split reads and filtered spanning reads)
  * **`*final_split_1.txt, final_split_2.txt`** (paired-end reads in fastq format: combine discordant and singlton split reads after filtering out unspecific read mapping at the genome level)
  * `final_spanning_1.txt, final_spanning_1.txt` (paired-end reads in fastq format: spanning reads after filtering out unspecific mapping at the genome level)
  * **`*final_spanning_noclip.sorted.bam`** (show spanning reads mapped to scaffold sequence, done by hisat2 no-splicing alignment model) # If there are no spanning reads, the "final_spanning_noclip.sorted.bam" is not present.
  * **`*final_split_noclip.sorted.bam`** (show filtered discordant and singlton split reads mapped to scaffold sequence, done by hisat2 no-splicing alignment model) # If there are no discordant/singlton split reads, the "final_split_noclip.sorted.bam" is not present.
  * `tmp` folder (if users want to look at more detail of processing steps):
  
    * `hisats_noclip.sam` (SAM format of "hisats_noclip.sorted.bam")
    * `spanning_sec.sam` (align "spanning_1.txt, spanning_2.txt" to genome reference, done by hisat2 splicing alignment model)
    * `discordant_split_sec.sam` (align "discordant_split_1.txt, discordant_split_2.txt" to genome reference, done by hisat2 splicing alignment model)
    * `singlton_split_sec.sam` (align "singlton_split_1.txt, singlton_split_2.txt" to genome reference, done by hisat2 splicing alignment model)
    
*: files with bold name are most important for users

## 4. Summarise the number of spanning and split reads across a cohort of samples (running evaluate.pl)
  4.1 Look at running parameters:
      perl evaluate.pl --help
      
  4.2 An example of running:
      Users have to create a directory that contains output folders (from `select_read.pl`) of the samples for summarizing. For instance in "examples" directory, `mkdir RCC1_ABHD12B_new && cp -r output RCC1_ABHD12B_new`, then run `evaluate.pl` as follows.
      
      perl evaluate.pl \
      
      --input ~/examples/RCC1_ABHD12B_new \
      # Set the input path
      
      --output ~/examples/RCC1_ABHD12B_summary \
      # Set the output directory of running "evaluate.pl"
  
  4.3 Summary of output directory run by `evaluate.pl`
  
     For example: ~/examples/RCC1_ABHD12B_summary
      
     * **`summary.txt`** (the number of read support [discordant_split, singlton_split; spanning] for breakpoint scaffold sequence across samples; statistics test for mapping bias to scaffold sequence [p<.05 indicates bias])
     
     * **`alt_1`** (concatenate split reads across all samples, and make alignment)
     
        |--- "All_sample_split_1.txt, All_sample_split_2.txt" (concatenate discordant/singlton split reads across all samples)
        |--- "All_sample_noclip.sorted.bam" (align concatenated split reads to breakpoint scaffold sequence, done by hisat2 no-splicing alignment model)
        |--- "p.value" (Fisher exact test for mapping bias to upstream/downstream of scaffold sequence)
        
## 5. Installation/Running via Docker

  5.1 [Install the Docker engine](https://docs.docker.com/engine/installation/) in your OS platform
  - installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/) 
  - installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/) 
  - installing [Docker on Windows](https://docs.docker.com/docker-for-windows/) (NOTE: We have not yet done enough testing on the Windows platform, so we would like to recieve more feedback on it)

  5.2 Allocate computational resource to docker, e.g.
  - Memory: min 4GB for ScaR running
  - CPUs: 4 (Users need to do setting based on their own hardwares)
  - Swap: 1GB (ScaR does not need a large memory for running, so keep a low Swap space)
 
  5.3 Pull / Build ScaR engine image

  5.3.1 Reqiurements
  - For Linux and Mac users, root privilege is needed. If you are non-root user, please refer to [this setting](https://docs.docker.com/install/linux/linux-postinstall/). 
  - Make sure that the internet is always accessible during Pulling / Building process (The implementation of dockerising is not suitable for TSD at this moment)
  
  5.3.2 Pull image from Docker Hub/Cloud repositories
  - Users can pull the ScaR engine image directly from DockerHub (approx 7.4Gb) which has been built and pushed to Docker Hub/Cloud repositories in advance. Run `docker pull senzhao/scar:latest`. After that, check the image by typing `docker images`
  
  5.3.2 Build image from docker container (optional)
  - If users would like to build the ScaR engine image instead of pulling it from Docker Hub, just download the soruce code and change to directory `cd ~/ScaR-master`, and then run `docker build --rm -t senzhao/scar:latest -f Dockerfile_ubunta .` (If the building process is not successful, please try another `docker build --rm -t senzhao/scar:latest -f Dockerfile_conda .`). NOTE: building is a long process (around 1-2 hours, dependent on network condition) and also needs a disk space with at leat free 50G.
  - After building is done, check the images by typing `docker images`
    
  5.4 Run ScaR engine image
  
  5.4.1 Usage - display all parameters: `docker run -t --rm senzhao/scar perl /ScaR/select_read.pl`
  
  5.4.2 Run an example using the data in the "examples" directory: 
  
  ```bash
  docker run -t --rm -v /input_data_path/examples:/data senzhao/scar perl /ScaR/select_read.pl --p 4 \
              --first input/raw_1.fastq \
              --second input/raw_2.fastq \
              --geneA RCC1 --geneB ABHD12B --trimm 0 \
              --scaffold input/RCC1_ABHD12B_scaff_seq.fa \
              --input /reference \
              --output output
  ```
  * "perl /ScaR/select_read.pl" - the running path of ScaR in docker image (keep it as this path and never make changes)
  * "/input_data_path/examples" - set the full path of input directory (it contains both raw reads and scaffold sequence files). Users need to define a new path for their own data.
  * "input/raw_1.fastq" - set the path of the R1 reads file in the input directory (relative path referring to input directory /input_data_path/examples). Users need to define file name of their own raw R1 reads. 
  * "input/raw_2.fastq" - set the path of the R2 reads file in the input directory (relative path referring to input directory /input_data_path/examples). Users need to define file name of their own raw R2 reads.
  * "input/RCC1_ABHD12B_scaff_seq.fa" - set the path of scaffold sequence file (relative path referring to input directory /input_data_path/examples). Users need to define file name of their own scaffold sequence. 
  * "/reference" - set the path of reference and annotation files (NOTE: keep it as "/reference" and never make changes)
  * "output" - set the output of ScaR running (relative path referring to the directory /input_data_path/examples).
  
  5.4.3 Run an example of `evaluate.pl` for summarizing the number of spanning and split reads across a cohort of samples: 
  Users have to create a directory that contains output folders (from `select_read.pl`) of the samples for summarizing. For instance in "examples" directory, `mkdir RCC1_ABHD12B_new && cp -r output RCC1_ABHD12B_new`, then run `evaluate.pl` as follows.
  
  ```bash
  docker run -t --rm -v /input_data_path/examples:/data senzhao/scar perl /ScaR/evaluate.pl \
      --input RCC1_ABHD12B_new --output RCC1_ABHD12B_summary
  ```

## 6. Reference
1. Kim D, Langmead B, and Salzberg SL, HISAT: a fast spliced aligner with low memory requirements. Nature Methods 12, 357-360 (2015). [DIO:10.1038/nmeth.3317](http://www.nature.com/nmeth/journal/v12/n4/full/nmeth.3317.html)
2. Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9 (2009). [DOI: 10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp352)

## Contact

t.cytotoxic@gmail.com

