ScaR package
========================

Use scaffold re-aligning approach to detect the prevalence and recurrence of known fusion transcripts across samples

![alt text]()

## 1. Requirements (before running program)
  
  1.1 Perl version >= 5.10.2
  
  1.2 HISAT2 version 2.1.0 (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip)
      
      The binary files have been integrated in ~/bin/hisat2-2.1.0/, please add the path to linux environment variables before running: 
        PATH=$PATH:/where_is_path/ScaR/bin/hisat2-2.1.0/
        export PATH

  1.3 HISAT2 aligner HGFM index (genome reference plus transcripts based on Ensembl GRCh38 version) 
      
      Users have to download the index files in ~/reference following the command:
        cd ~/reference
        wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz"
        tar -vxf grch38_tran.tar.gz
        mv grch38_tran/genome_tran.* .
  
  1.4 Samtools version >= 1.3 (https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2)
      
      If samtools has not been installed in the system, users have to download and install it locally.
      Then, add the path to linux environment variables before running:
        PATH=$PATH:/where_is_path_samtools
        export PATH
 
  1.5 R version >= 3.0.3 (https://cran.r-project.org)
      
      If R has not been installed in the system, users have to donwload and install it locally.
      Then, add the' path to linux environment variables before running:
        PATH=$PATH:/where_is_path/R
        export PATH
    
  1.6 Genomic data and annotations
  
      Type the command and download these files
        cd ~/reference
        wget "http://folk.uio.no/senz/Gene_hg38.txt" -- gene annotation file
        wget "http://folk.uio.no/senz/ensembl_transcript.fa" -- the human transcriptome sequences annotated from ensembl database (Ensembl Archive Release 89)
        wget "http://folk.uio.no/senz/gencode_transcript.fa" -- the human transcriptome sequences annotated from GENCODE database (Release version 27)
        wget "http://folk.uio.no/senz/ucsc_transcript.fa" -- the human transcriptome sequences annotated from UCSC database (Release date: Jan 2018)
        wget "http://folk.uio.no/senz/ucsc_refGene.txt"
        
  1.7 Set the path of Perl libraries to environment variables
  
      PERL5LIB="$PERL5LIB:/where_is_path/ScaR/lib"
      export PERL5LIB
      NOTE: we recommend that users add the path of perl libraries to .bashrc, then "source .bashrc"
        
## 2. Usage of the "select_read.pl"
  2.1 See running parameters
      
      perl select_read.pl --help

  2.2 An example of running:
      
      perl select_read.pl \
      
      --first ~/examples/input/raw_1.fastq \ 
      # Raw or compressed fastq (.fastq.gz) file for the first end of paired-end reads
      
      --second ~/examples/input/raw_2.fastq \
      # Raw or compressed fastq (.fastq.gz) file for the second end of paired-end reads
      
      --geneA RCC1 --geneB ABHD12B \
      # Fusion partner gene names (Refseq gene symbol and Ensembl id are accepted currently)
      
      --anchor 6 \ 
      # (default: 6)
      # Set the length of anchor. The minimum number of bases is required to match to geneA/geneB region in the scaffold sequence.
      
      --trimm 0 \
      # (default: 0)
      # Set whether the input fastq reads are trimmed (1) or not (0)
      
      --length 48 \
      # Set the maximum value of fastq read length, this option is only active when raw fastq reads are trimmed (--trimm 1)
      
      --trans_ref ensembl
      # (default: ensembl)
      # The setting of annotation resources, users could choose other options (e.g. "gencode" or "ucsc").
      
      --p 8 \ 
      # (default: 8)
      # The number of threads for running in parallel. 
      
      --anno ~/reference/ \
      # Set the path of genomic data and annotation
      
      --output ~/examples/output/ \
      # Output directory
      
      --scaffold ~/examples/input/RCC1_ABHD12B_scaff_seq.fa \
      # Fusion scaffold sequences, users can extract the breakpoint sequences from outputs of de novo fusion finders (e.g. deFuse, fusioncatcher, SOAPfuse). Candidate sequences in fasta format are accepted. 
      # For instance:
      #	>alt_0
      #	XXXXXXXXXXXXXXXXXX|YYYYYYYYYYYYYYYY
      #	>alt_1
      #	XXXXXXXXXXXXXXXXXX*YYYYYYYYYYYYYYY
      # NOTE: 1. It only accepts '*' or '|' as a separator for breakpoint sequences. 'XXXXXXXXXXX' and 'YYYYYYYYY' correspond to the sequences from geneA and geneB, respectively.
      #       2. In order to ensure the specificity of breakpoint sequences matching to the reference, we recommend that 'XXXXXXXX' and 'YYYYYYYY' should be at least 20 bp.
      #       3. In general, the breakpoint sequences are composed of cDNAs (i.e. exon region). If users would like to detect the fusion sequences including intron/intergenic region, they have to set user-defined reference sequences , please see the usage of parameter "--user_ref".

      --user_ref ~/upstream.fasta
      # User-defined reference sequences, user can specify transcript reference sequences (or genomic sequences) in fasta format which are not present in the default database.
      # For instance:
      # >RP11-599B13.3|alternative1
      # CTTTGTGTCTTTGTCTTTATTTCTTTTCTCATTCCCTCGTCTCCACCGGGAAGGGGAGAGCCTGCGGGTGGTGTATCAGGCAGGTTCCCCTACATCTTTGGCACCCAACAC
      # NOTE: 'RP11-599B13.3' is the gene name and should be identical to the input of gene partner names (either GeneA or GeneB); 'alternative1' is the transcript name (please avoid using symbol '_' in user-defined transcript name). Make sure both 'RP11-599B13.3' and 'alternative1' are present together, and separated by '|'.
    
      
## 3. Output results
  For example: ~/examples/output/
  
  * `scaffold_*_seq.fa` (cDNA sequences of geneA and geneB, and breakpoint sequence of scaffold in fasta format)
  * `hisats_noclip.sorted.bam` (BAM file: done by no-splicing alignment model)
  * `discordant_split_1.txt, discordant_split_2.txt` (paired-end reads in fastq format: one-end maps to scaffold; the other maps to cDNA sequences of geneA/geneB)
  * `singlton_split_1.txt, singlton_split_2.txt` (paired-end reads in fastq format: one-end maps to scaffold; the other has no mapping to cDNA sequences of geneA/geneB)
  * `spanning_1.txt, spanning_2.txt` (paired-end reads in fastq format: one-end maps to cDNA sequence of geneA; the other maps to cDNA sequence of geneB)
  * `read_mapped_info` (mapping summary of discordant/singlton split reads and spanning reads)
  * **`*final_read_mapped_info`** (mapping summary of filtered discordant/singlton split reads and filtered spanning reads)
  * **`*final_split_1.txt, final_split_2.txt`** (paired-end reads in fastq format: merge discordant and singlton split reads after filtering out unspecific mapping)
  * **`*final_spanning_1.txt, final_spanning_1.txt`** (paired-end reads in fastq format: spanning reads after filtering out unspecific mapping)
  * **`*final_spanning_noclip.sorted.bam`** (filtered spanning reads mapped to scaffold sequence). If there are no spanning reads, the "final_spanning_noclip.sorted.bam" is not present.
  * **`*final_split_noclip.sorted.bam`** (filtered discordant and singlton split reads mapped to scaffold sequence). If there are no discordant/singlton split reads, the "final_split_noclip.sorted.bam" is not present.
  * **`*summary_read_mapping_support.txt`** (comparison of read supports for discordant split and spanning reads before/after filtering out unspecific mapping)
    
*: files with bold name are important for users

## 4. Summarise the number of spanning and split reads across a cohort of samples (running evaluate.pl)
  4.1 See running parameters:
      perl evaluate.pl --help
      
  4.2 An example of running:
      Users have to create a directory that contains outputs from `select_read.pl` for summarizing. For instance in "examples" directory, if RCC1_ABHD12B_new folder is not present, `mkdir RCC1_ABHD12B_new && cp -r output RCC1_ABHD12B_new/`, then run `evaluate.pl` as follows.
      
      perl evaluate.pl \
      
      --input ~/examples/RCC1_ABHD12B_new \
      # Set the input path
      
      --output ~/examples/RCC1_ABHD12B_summary \
      # Set the output directory of running "evaluate.pl"
  
  4.3 Summary of output directory run by `evaluate.pl`
  
     For example: ~/examples/RCC1_ABHD12B_summary
      
     * **`summary.txt`** (the number of read support [discordant_split, singlton_split; spanning] for breakpoint scaffold sequence across a cohort of samples; statistics test for mapping distribution bias to scaffold sequence [p<.05 indicates bias])
     
     * **`alt_0`** (concatenate split reads across a cohort of samples, and make an alignment)
     
        |--- "All_sample_split_1.txt, All_sample_split_2.txt" (paired-end discordant/singlton split reads concatenated across all samples)
        |--- "All_sample_noclip.sorted.bam" (align the concatenated split reads to scaffold sequence, user can visualize this bam file uisng scaffold_ENST00000373833_49_ENST00000337334_48_seq.fa as reference by IGV)
        |--- "p.value" (fisher exact test for mapping distribution bias to upstream/downstream of scaffold sequence)
        
## 5. Installation/Running via Docker

  5.1 [Install the Docker engine](https://docs.docker.com/engine/installation/) in your OS platform
  - installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/) 
  - installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/) 
  - installing [Docker on Windows](https://docs.docker.com/docker-for-windows/) (NOTE: We have not yet done enough testing on the Windows platform, so we would like to recieve more feedback on it)

  5.2 Allocate computational resource to docker, e.g.
  - Memory: min 4GB for ScaR running
  - CPUs: 4 (Users have to set up based on their own hardwares)
  - Swap: 1GB (ScaR does not need a large memory for running, so keep a small mount of Swap space)
 
  5.3 Pull / Build ScaR engine image

  5.3.1 Reqiurements
  - For Linux and Mac users, root privilege is needed. If you are non-root user, please refer to [this setting](https://docs.docker.com/install/linux/linux-postinstall/). 
  - Make sure that the internet is always accessible during Pulling / Building process.
  
  5.3.2 Pull image from Docker Hub/Cloud repositories
  - Users can pull the ScaR engine image directly from Docker Hub (approx 7.4Gb) which has been built and pushed to Docker Hub/Cloud repositories in advance. Run `docker pull senzhao/scar:latest`. After that, check the image by typing `docker images`
  
  5.3.2 Build image from docker container (optional)
  - If users would like to build the ScaR engine image instead of pulling it from Docker Hub, just download the soruce code and change to directory `cd ~/ScaR-master`, and then run `docker build --rm -t senzhao/scar:latest -f Dockerfile_ubunta .` (If this building process is not successful, please try another `docker build --rm -t senzhao/scar:latest -f Dockerfile_conda .`). NOTE: building is a long process (around 1-2 hours, dependent on network condition) and also needs a disk space with at least free 50G.
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
              --anno /reference \
              --output output
  ```
  * "perl /ScaR/select_read.pl" - the running path of ScaR in docker image (NOTE: keep it as this path)
  * "/input_data_path/examples" - set the full path of input directory (it contains both raw reads and scaffold sequence files). Users need to define a new path for their own data.
  * "input/raw_1.fastq" - set the path of the R1 reads file in the input directory (relative path referring to input directory /input_data_path/examples). users need to define file name of their own R1 reads. 
  * "input/raw_2.fastq" - set the path of the R2 reads file in the input directory (relative path referring to input directory /input_data_path/examples). users need to define file name of their own R2 reads.
  * "input/RCC1_ABHD12B_scaff_seq.fa" - set the path of scaffold sequence file (relative path referring to input directory /input_data_path/examples). users need to define file name of their own scaffold sequences. 
  * "/reference" - the path of reference and annotation files in docker image (NOTE: keep it as "/reference")
  * "output" - set the output of ScaR running (relative path referring to the directory /input_data_path/examples).
  
  5.4.3 Run an example of `evaluate.pl` for summarizing the number of spanning and split reads across a cohort of samples: 
  Users have to create a directory that contains output from `select_read.pl` for summarizing. For instance in "examples" directory, if RCC1_ABHD12B_new folder is not present, `mkdir RCC1_ABHD12B_new && cp -r output RCC1_ABHD12B_new/`, then run `evaluate.pl` as follows.
  
  ```bash
  docker run -t --rm -v /input_data_path/examples:/data senzhao/scar perl /ScaR/evaluate.pl \
      --input RCC1_ABHD12B_new --output RCC1_ABHD12B_summary
  ```

## 6. Reference
1. Kim D, Langmead B, and Salzberg SL, HISAT: a fast spliced aligner with low memory requirements. Nature Methods 12, 357-360 (2015). [DIO:10.1038/nmeth.3317](http://www.nature.com/nmeth/journal/v12/n4/full/nmeth.3317.html)
2. Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9 (2009). [DOI: 10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp352)
3. Zhao S*., Andreas M Hoff*, Rolf I Skotheim, ScaR - A tool for sensitive detection of known fusion transcripts: Establishing prevalence of fusions in testicular germ cell tumors. BioRix [https://www.biorxiv.org/content/10.1101/518316v1](https://www.biorxiv.org/content/10.1101/518316v1) (* equal contribution)

## Contact

t.cytotoxic@gmail.com

