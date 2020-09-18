ScaR package
========================

Use scaffold re-aligning approach to detect the prevalence and recurrence of known fusion transcripts across samples

![My image](https://github.com/senzhaocode/ScaR/blob/master/examples/demo.png?raw=true)

## 1. Requirements (before running program)
  
  1.1 Perl version >= 5.10.2
  
  1.2 HiSAT2 v2.1.0 (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip)
      
      The binary executable files have been integrated in ~/bin/hisat2-2.1.0/, please add a working path to Linux environment variables before running: 
        PATH=$PATH:/where_is_path/ScaR/bin/hisat2-2.1.0/
        export PATH

  1.3 HiSAT2 aligner HGFM index files (genome reference plus transcripts based on Ensembl GRCh38 version) 
      
      Users have to download the index files in ~/reference following the settings:
        cd ~/reference
        wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz"
        tar -vxf grch38_tran.tar.gz
        mv grch38_tran/genome_tran.* .
  
  1.4 Samtools version >= 1.3 (https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2)
      
      If samtools is not available in the system, users have to download and install it locally. 
        PATH=$PATH:/where_is_path_samtools
        export PATH
 
  1.5 R version >= 3.0.3 (https://cran.r-project.org)
    
  1.6 Genomic data and annotations
  
        cd ~/reference
        wget "http://folk.uio.no/senz/GRCh38.primary_assembly.genome.fa" #-- whole genome sequences (build GRCh38 version)
        wget "http://folk.uio.no/senz/Gene_hg38.txt" #-- gene annotation file
        wget "http://folk.uio.no/senz/ensembl_transcript.fa" #-- the human transcriptome sequences annotated from ensembl database (Ensembl Archive Release 89)
        wget "http://folk.uio.no/senz/gencode_transcript.fa" #-- the human transcriptome sequences annotated from GENCODE database (Release version 27)
        wget "http://folk.uio.no/senz/ucsc_transcript.fa" #-- the human transcriptome sequences annotated from UCSC database (Release date: Nov 2018)
        wget "http://folk.uio.no/senz/ucsc_refGene.txt"
        
  1.7 Set Perl libraries to environment variables
  
      PERL5LIB="$PERL5LIB:/where_is_path/ScaR/lib"
      export PERL5LIB
      NOTE: we recommend that users add the path of perl libraries to .bashrc, then "source .bashrc"

## 2. Optional dependencies

  If users would like to use an optional aligner STAR instead of the default HiSAT2:
  
  2.1 STAR v2.7.2d (https://github.com/alexdobin/STAR/archive/2.7.2d.tar.gz)
  
      The binary executable files have been integrated in ~/bin/STAR-2.7.2d/, please add a working path to Linux environment variables before running:
        PATH=$PATH:/where_is_path/ScaR/bin/STAR-2.7.2d/
        export PATH
        
  2.2 STAR aligner SA index files (genome reference plus transcript annotation in gtf format)
  
      Users can download index files that were pre-build on basis of GRCh38 genome reference and Ensembl v89 transcript annotations)
        cd ~/reference
        wget "http://folk.uio.no/senz/STAR_index.tar.gz"
        tar -vxf STAR_index.tar.gz
        mv STAR_index/* .
        
        NOTE: Users can generate their own genome indexes with most latest assemblies and annotations.

## 3. Usage of the "select_read.pl"
  3.1 See running parameters
      
      perl select_read.pl --help

  3.2 An example of running:
      
      perl select_read.pl
      
      --first ~/examples/input/raw_1.fastq 
      # Raw *.fastq or compressed *.fastq.gz file for the 1st end of paired-end reads
      
      --second ~/examples/input/raw_2.fastq 
      # Raw *.fastq or compressed *.fastq.gz file for the 2nd end of paired-end reads
      
      --geneA RCC1 --geneB ABHD12B 
      # Fusion partner gene names (Refseq gene symbol and Ensembl id are accepted)
      
      --anchor 6 
      # (default: 6)
      # Set the length of anchor. The minimum number of bases is required to match to geneA/geneB region in the scaffold sequence.
      
      --trimm 0 
      # (default: 0)
      # Set whether the input fastq reads are trimmed (1) or not (0)
      
      --length 48 
      # Set the maximum length of reads, this option is only available when raw fastq reads are trimmed (--trimm 1)
      
      --transAlign
      # (default: hisat2)
      # Set an aligner to map reads at transcriptome level using no-splicing mode (Options: hisat2 or star)
      
      --genomeAlign
      # (default: hisat2)
      # Set an aligner to map reads at genome level using splicing mode (Options: hisat2 or star)
      
      --trans_ref ensembl
      # (default: ensembl)
      # Set the resource of transcript sequence data (Options: gencode, ensembl or ucsc).
      
      --p 8 
      # (default: 8)
      # The number of threads for running in parallel. 
      
      --anno ~/reference/ 
      # Set the directory path of genomic, transcriptomic sequences and annotations
      
      --output ~/examples/output/ 
      # Output directory
      
      --scaffold ~/examples/input/RCC1_ABHD12B_scaff_seq.fa 
      # A list of fusion scaffold sequences in fasta format (if --scaffold is active, --coordinate has to be inactivated)
      # For instance:
      #	>alt_0
      #	XXXXXXXXXXXXXXXXXX|YYYYYYYYYYYYYYYY
      #	>alt_1
      #	XXXXXXXXXXXXXXXXXX*YYYYYYYYYYYYYYY
      # NOTE: 1. '*' or '|' is accepted as a separator for breakpoint sequences.
      #       2. To ensure the specificity of breakpoint sequences matching to the reference, we recommend that the lengths of 'XXXXXXXX' and 'YYYYYYYY' have to be at least 20 bp.
      #       3. In general, the breakpoint sequences are composed of cDNAs (i.e. exon region). If users would like to detect the fusion sequences including intron/intergenic region, they have to set user-defined reference sequences, please see the usage of parameter "--user_ref".

      --coordinate "chr1:34114119|chr2:65341523,chr1:3412125|chr2:65339145" or "chr1:34114119:+|chr2:65341523:-,chr1:3412125|chr2:65339145"
      # Set genomic junction coordinates (build GRCh38) of breakpoint sites and strand directions (optional input) for GeneA and GeneB, e.g. chr1:34114119:+ and chr2:65341523:- correpsond to the chromosome names, genomic breakpoints and strand directions (optional) of GeneA and GeneB, respectively. (if --coordinate is active, --scaffold has to be inactivated)
      
      --user_ref ~/upstream.fasta
      # User-defined reference sequences, users can specify transcript reference sequences (or genomic sequences) in fasta format which are not present in the provided annotation databases.
      # For instance:
      # >RP11-599B13.3|alternative1
      # CTTTGTGTCTTTGTCTTTATTTCTTTTCTCATTCCCTCGTCTCCACCGGGAAGGGGAGAGCCTGCGGGTGGTGTATCAGGCAGGTTCCCCTACATCTTTGGCACCCAACAC
      # NOTE: 'RP11-599B13.3' is the gene name and has to be identical to the input of gene partner names; 'alternative1' is the transcript name (please avoid using the symbol '_ $ % & # * @ ^ ? + ! < > | / \' in user-defined transcript name). Make sure both 'RP11-599B13.3' and 'alternative1' are together, and separated by '|'.
    
      
## 4. Output results
  For example: ~/examples/output/
  
  * `scaffold_*_seq.fa` (cDNA sequences of geneA and geneB, and breakpoint sequence of scaffold in fasta format)
  * `discordant_split_1.txt, discordant_split_2.txt` (paired-end reads in fastq format: one-end maps to scaffold; the other maps to geneA/geneB)
  * `singlton_split_1.txt, singlton_split_2.txt` (paired-end reads in fastq format: one-end maps to scaffold; the other has no mapping to geneA/geneB)
  * `spanning_1.txt, spanning_2.txt` (paired-end reads in fastq format: one-end maps to geneA; the other maps to geneB)
  * `read_mapped_info` (mapping summary of discordant/singlton split reads and spanning reads)
  * **`*final_read_mapped_info`** (mapping summary of discordant/singlton split and spanning reads after filtering out unspecific mapping)
  * **`*final_split_1.txt, final_split_2.txt`** (merge discordant and singlton split reads in fastq format after filtering out unspecific mapping)
  * **`*final_spanning_1.txt, final_spanning_1.txt`** (spanning reads in fastq format after filtering out unspecific mapping)
  * **`*final_spanning_noclip.sorted.bam`** (filtered spanning reads mapped to GeneA/GeneB/scaffold sequence). If there are no spanning reads, "final_spanning_noclip.sorted.bam" is absent.
  * **`*final_split_noclip.sorted.bam`** (filtered discordant/singlton split reads mapped to GeneA/GeneB/scaffold sequence). If there are no discordant/singlton split reads, "final_split_noclip.sorted.bam" is absent.
  * **`*summary_read_mapping_support.txt`** (the number and fraction of discordant split and spanning reads that show a unique alignment to GeneA, GeneB and scaffold, after filtering out unspecific mapped reads)
    
*: files with bold name are important for users

## 5. Summarize the number of spanning and split reads across a cohort of samples (run evaluate.pl)
  5.1 See running parameters:
      perl evaluate.pl --help
      
  5.2 An example of running:
      Users have to make a directory that contains outputs from `select_read.pl` for summarizing. For instance in "examples" directory, if RCC1_ABHD12B_new folder is not present, `mkdir RCC1_ABHD12B_new && cp -r output RCC1_ABHD12B_new/`, then run `evaluate.pl` as follows:
      
      perl evaluate.pl 
      
      --input ~/examples/RCC1_ABHD12B_new 
      # Set the input path
      
      --output ~/examples/RCC1_ABHD12B_summary 
      # Set the output directory of running "evaluate.pl"
  
  5.3 Summarize output directory of `evaluate.pl`
  
     For example: ~/examples/RCC1_ABHD12B_summary
      
     * **`summary.txt`** (the number of read support [discordant_split, singlton_split; spanning] for breakpoint scaffold sequence across a cohort of samples; statistics test for mapping distribution bias to scaffold sequence)
     
     * **`alt_0`** (concatenate all split reads across a cohort of samples, and align to scaffold sequence)
     
        |--- "All_sample_split_1.txt, All_sample_split_2.txt" (paired-end discordant/singlton split reads concatenated across all samples)
        |--- "All_sample_noclip.sorted.bam" (align concatenated split reads to scaffold sequence, and user can visualize this bam file uisng scaffold_ENST00000373833_49_ENST00000337334_48_seq.fa as reference by IGV)
        |--- "p.value" (fisher exact test and chi-squred test for mapping distribution bias to upstream/downstream of scaffold sequence)
        
## 6. Installation/Running via Docker

  6.1 [Install the Docker engine](https://docs.docker.com/engine/installation/) in your OS platform
  - installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/) 
  - installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/) 
  - installing [Docker on Windows](https://docs.docker.com/docker-for-windows/) (NOTE: We have not yet done enough testing on the Windows platform, so we would like to recieve more feedback on it)

  6.2 Allocate computational resource to docker, e.g.
  - Memory: min 4GB for ScaR running
  - CPUs: 4 (Users have to set up based on their own hardwares)
  - Swap: 1GB (ScaR does not need a large memory for running, so keep a small mount of Swap space)
 
  6.3 Pull / Build ScaR engine image

  6.3.1 Requirements
  - For Linux and Mac users, root privilege is needed. If you are non-root user, please refer to [this setting](https://docs.docker.com/install/linux/linux-postinstall/). 
  
  6.3.2 Pull image from Docker Hub/Cloud repositories
  - Users can pull the ScaR engine image directly from Docker Hub (approx 9Gb) which has been built and pushed to Docker Hub/Cloud repositories in advance. Run `docker pull senzhao/scar:latest`. After that, check the image by typing `docker images`
  
  6.3.3 Build image from docker container (optional)
  - If users would like to build the ScaR engine image instead of pulling it from Docker Hub, just download the soruce code and change to directory `cd ~/ScaR-master`, and then run `docker build --rm -t senzhao/scar:latest -f Dockerfile_ubunta .`. NOTE: building may take a long process (around 30 mins) and also needs a disk image size with at least 50G.
  - After building is done, check the images by typing `docker images`
    
  6.4 Run ScaR engine image
  
  6.4.1 Usage - display all parameters: `docker run -t --rm senzhao/scar perl /ScaR/select_read.pl`
  
  6.4.2 Run an example using the data in the "examples" directory: 
  
  ```bash
  docker run -t --rm -v /input_data_path/examples:/data senzhao/scar perl /ScaR/select_read.pl --p 4 \
              --first input/raw_1.fastq \
              --second input/raw_2.fastq \
              --geneA RCC1 --geneB ABHD12B --trimm 0 \
              --scaffold input/RCC1_ABHD12B_scaff_seq.fa \
              --anno /reference \
              --output output
              
  docker run -t --rm -v /input_data_path/examples:/data senzhao/scar perl /ScaR/select_read.pl --p 4 \
              --first input/raw_1.fastq \
              --second input/raw_2.fastq \
              --geneA RCC1 --geneB ABHD12B --trimm 0 \
              --coordinate "chr1:28508159|chr14:50901829" \
              --anno /reference \
              --output output
  ```
  
  ***Some notes for running ScaR using STAR aligner under Docker framework:***
  
  - As STAR genome index files are very large, we do not pack them into container and create a huge image (it can be harder to distribute and use). The best way to load the STAR genome index files in host machine is to use a "bind mount" solution to attach the volume to the container. For example, if the storage directory in host machine is `/input_data_path/examples`, 
    
        cd /input_data_path/examples
        wget "http://folk.uio.no/senz/STAR_index.tar.gz"
        tar -vxf STAR_index.tar.gz
    
  - If user would like to use the index files generated by theirselves, please make a new folder named `STAR_index` (always keep this name) and transfer all index files (e.g. SA, SAindex, Genome) within `/input_data_path/examples/STAR_index`.

```bash
  docker run -t --rm -v /input_data_path/examples:/data senzhao/scar perl /ScaR/select_read.pl --p 4 \
              --first input/raw_1.fastq \
              --second input/raw_2.fastq \
              --geneA RCC1 --geneB ABHD12B --trimm 0 \
              --transAlign star \
              --genomeAlign star \
              --scaffold input/RCC1_ABHD12B_scaff_seq.fa \
              --anno /reference \
              --output output
``` 
    
  * "perl /ScaR/select_read.pl" - the running path of ScaR in docker image (NOTE: keep it as "/ScaR/select_read.pl").
  * "/input_data_path/examples" - set the full path of input directory (it contains both raw reads and scaffold sequence files) in host machine.
  * "input/raw_1.fastq" - set the path of R1-end read file (relative path to input directory /input_data_path/examples in host machine).
  * "input/raw_2.fastq" - set the path of R2-end read file (relative path to input directory /input_data_path/examples in host machine).
  * "input/RCC1_ABHD12B_scaff_seq.fa" - set the path of scaffold sequence file (relative path to input directory /input_data_path/examples in host machine). 
  * "/reference" - the path of reference and annotation files in docker image (NOTE: keep it as "/reference")
  * "output" - set the output of ScaR running (relative path to the directory /input_data_path/examples in host machine).
  
  6.4.3 Run an example of `evaluate.pl` for summarizing the number of spanning and split reads across a cohort of samples: 
  Users have to make a directory that contains outputs from `select_read.pl` for summarizing. For instance in "examples" directory, if RCC1_ABHD12B_new folder is not present, `mkdir RCC1_ABHD12B_new && cp -r output RCC1_ABHD12B_new/`, then run `evaluate.pl` as follows.
  
  ```bash
  docker run -t --rm -v /input_data_path/examples:/data senzhao/scar perl /ScaR/evaluate.pl \
      --input RCC1_ABHD12B_new --output RCC1_ABHD12B_summary
  ```

## 7. Reference
1. Kim D, Langmead B, and Salzberg SL, HISAT: a fast spliced aligner with low memory requirements. Nature Methods 12, 357-360 (2015). [DIO:10.1038/nmeth.3317](http://www.nature.com/nmeth/journal/v12/n4/full/nmeth.3317.html)
2. Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9 (2009). [DOI: 10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btp352)
3. Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15–21. [doi:10.1093/bioinformatics/bts635](https://academic.oup.com/bioinformatics/article/29/1/15/272537)
4. Zhao S*., Andreas M Hoff*, Rolf I Skotheim, ScaR—a tool for sensitive detection of known fusion transcripts: establishing prevalence of fusions in testicular germ cell tumors. NAR Genomics and Bioinformatics, Volume 2, Issue 1, March 2020, lqz025, [https://doi.org/10.1093/nargab/lqz025](https://academic.oup.com/nargab/article/2/1/lqz025/5701460) (* equal contribution)

## Contact

t.cytotoxic@gmail.com

