#!/usr/bin/perl -w
###########################################################################################################################
# Select_read.pl is perl script to detect the recurrence of known fusion transcripts using scaffold realignment strategy 
# Before you start it, please read manual first: perl select_read.pl --help
# The script was created in Nov. 11 2016, 
# The second modification in Nov 26 2017 (update annotation), 
# The third modification in Dec 2nd, 2017 (resolve anchor problem for read mapping to scafflod),
# The forth modification in June, 2018 (add --user_ref parameter, allow user-defined transcript sequences as input),
###########################################################################################################################
use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Extract;
use First_output;
use Fusion_gene;
use Genome_align;
use Check;

	# setup option parameters in the script command line -- for RNA-seq analyses
	my @usage;
	push @usage, "Usage: ".basename($0)." [options]\n";
	push @usage, "Retrieve discordant (and singlton) split reads that support a given fusion breakpoint.\n";
	push @usage, "	--help		Displays this information\n";
	push @usage, "	--first		Raw fastq file or compressed fastq (.fastq.gz) file for 1st end of paired-end reads\n";
	push @usage, "	--second	Raw fastq file or compressed fastq (.fastq.gz) file for 2nd end of paired-end reads\n";
	push @usage, "	--geneA		Name of upstream gene partner (Gene_symbol or Ensembl_id is accpeted)\n";
	push @usage, "	--geneB		Name of downstream gene partner (Gene_symbol or Ensembl_id is accepted)\n";
	push @usage, "	--scaffold	A list of fusion scaffold sequences in fasta format, e.g\n", 
		     "			>scaffold1\n",
		     "			GCTCTATGAAATTGCA|AACAAAGAGAGGGTCA\n",
		     "			>scaffold2\n",
		     "			GCTCTATGAAATTGCA*AACAAAGAGAGGGTCA\n";
	push @usage, "	--anno		Set the directory of annotation files: e.g. /script_path/data/\n";
	push @usage, "	--output	Output Directory\n";
	push @usage, "	--anchor	The length of anchor for read mapping to breakpoint (defalut: 6)\n";
	push @usage, "	--trimm		Set whether input fastq reads are trimmed (1) or not (0, as default value)\n";
	push @usage, "	--length	Set the maximum length of fastq reads, this option is only available when raw fastq reads are trimmed (--trimm 1)\n";
	push @usage, "	--trans_ref	Set the resource of transcript sequence data (e.g. gencode, ensembl or ucsc, default: ensembl)\n";
	push @usage, "	--user_ref	Set the resource of user-defined transcript sequences (optional)\n";
	push @usage, "	--p		The number of threads, and make sure that it should be exact as the number of CPUs allocated in jobscript (defalut: 8)\n";
	my $help;
	my $fastq_1;
	my $fastq_2;
	my $geneA;
	my $geneB;
	my $scaffold;
	my $input;
	my $output;
	my $anchor;
	my $trimm;
	my $read_length;
	my $trans_ref;
	my $user_ref;
	my $cpus;

	GetOptions
	(
        	'help'        => \$help,
        	'first=s'     => \$fastq_1,
        	'second=s'    => \$fastq_2,
        	'geneA=s'     => \$geneA,
        	'geneB=s'     => \$geneB,
		'scaffold=s'  => \$scaffold,
		'anno=s'     => \$input,
        	'output=s'    => \$output,
        	'anchor=i'    => \$anchor,
		'trimm=i'     => \$trimm,
	 	'length=i'    => \$read_length,
	 	'trans_ref=s' => \$trans_ref,
		'user_ref=s'  => \$user_ref,
        	'p=i'         => \$cpus
	);
	not defined $help or die @usage;
	defined $fastq_1 or die @usage; if (! -e $fastq_1 ) { print "\n$fastq_1 is wrong path, please set valid path\n\n"; exit; }
	defined $fastq_2 or die @usage; if (! -e $fastq_2 ) { print "\n$fastq_2 is wrong path, please set valid path\n\n"; exit; }
	defined $scaffold or die @usage; if (! -e $scaffold ) { print "\n$scaffold is wrong path, please set valid path\n\n"; exit; }
	defined $geneA or die @usage;
	defined $geneB or die @usage;
	defined $output or die @usage;
	defined $input or die @usage; if (! -e "$input/ensembl_transcript.fa" ) { print "cDNA sequence file $input/ensembl_transcript.fa does not exist\n\n"; exit; } 
				if (! -e "$input/gencode_transcript.fa" ) { print "cDNA sequence file $input/gencode_transcript.fa does not exist\n\n"; exit; }
				if (! -e "$input/ucsc_transcript.fa" ) { print "cDNA sequence file $input/ucsc_transcript.fa does not exist\n\n"; exit; }
				if (! -e "$input/ucsc_refGene.txt" ) { print "Gene annotation file $input/ucsc_refGene.txt does not exist\n\n"; exit; }
				if (! -e "$input/Gene_hg38.txt" ) { print "Gene annotation file $input/Gene_hg38.txt does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.1.ht2" ) { print "\nHisat2 index file $input/genome_tran.1.ht2 does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.2.ht2" ) { print "\nHisat2 index file $input/genome_tran.2.ht2 does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.3.ht2" ) { print "\nHisat2 index file $input/genome_tran.3.ht2 does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.4.ht2" ) { print "\nHisat2 index file $input/genome_tran.4.ht2 does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.5.ht2" ) { print "\nHisat2 index file $input/genome_tran.5.ht2 does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.6.ht2" ) { print "\nHisat2 index file $input/genome_tran.6.ht2 does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.7.ht2" ) { print "\nHisat2 index file $input/genome_tran.7.ht2 does not exist\n\n"; exit; }
				if (! -e "$input/genome_tran.8.ht2" ) { print "\nHisat2 index file $input/genome_tran.8.ht2 does not exist\n\n"; exit; }
	if (! defined($anchor) ) { $anchor = 6; } # set min size to match margin of upstream / downstream in breakpoint sequence (default: 6 bp)
	if (! defined($read_length) ) { $read_length = 0; }
	if (! defined($trimm) ) { # judge the setting of $trimm parameter 
		$trimm = 0; # no trimming process
		$read_length = 0; # set initial value of read_lenght as 0
	} else {
		if ( $trimm == 0 ) {
			$read_length = 0;
		} else {
			if ( $read_length == 0 ) {
				print "The input read length is missing or invalid, pleaes input the maximum read length\n"; exit;
			}
		}
	}
	
	if (! defined($trans_ref) ) { # judge the setting of $trans_ref parameter
		$trans_ref = "ensembl"; 
	} else {
		if ( $trans_ref eq 'ensembl' ) {
		} elsif ( $trans_ref eq 'gencode' ) {
		} elsif ( $trans_ref eq 'ucsc' ) {
			if ( $geneA =~/^ENSG/ ) { print "Please use the gene symbol name as input when UCSC transcript sequences are implemented\n"; exit; }
			if ( $geneB =~/^ENSG/ ) { print "Please use the gene symbol name as input when UCSC transcript sequences are implemented\n"; exit; }
		} else {
			print "The transcript sequeuce resource does not exists, please use the valid resource (e.g. ensembl, gencode or ucsc)\n"; exit;
		}
	}

	if ( defined($user_ref) ) { # judge the setting of user defined transcript reference 
		if (! -e "$user_ref" ) {
			print "User defined transcript reference sequences do not exist, please set the valid path\n\n"; exit;
		}
	} else {
		$user_ref = ""; # no user-defined reference sequence presence
	}
	
	if (! defined($cpus) ) { $cpus = 8; } # set the default number of threads to run script (default: 8)

	# setup input environment variables
	my $genome_index = "$input/genome_tran"; # set Hisat2 index
	my $gene = "$input/Gene_hg38.txt"; # gene annotation files

	# check output folder status
	if ( -e $output ) { # if output folder present, please clean up first.
		if (! -z $output ) {
			`rm -r $output/*`; print "\nWARNINGS: $output exists, and the old files will be re-written\n\n";
		}
	} else {
		`mkdir $output`;
	}

print "############################################################################################\n";
print "# Step 0: quality control of  the input scaffold sequence / gene symbol name / read length #\n";
print "############################################################################################\n\n";
	# calculate the read length
	# Note: we recommend not trimming raw reads
	if (! $read_length ) { $read_length = Check::read_length($fastq_1, $fastq_2); }
	if (! $read_length ) { print "\nStep 0-1: Read length of fastq file is error, please set valid RNA-seq fastq file\n"; exit; }

	# judge the header name of fastq -- most likely three possibilities: " ", "/" and "end"
	my $header_sep = Check::judge_header($fastq_1);
	if ( $header_sep eq "end" ) { print "\nStep 0-2: The header of fastq_1 file is error, please make sure the fastq format is valid\n"; exit; }
	$header_sep = Check::judge_header($fastq_2);
	if ( $header_sep eq "end" ) { print "\nStep 0-2: The header of fastq_2 file is error, please make sure the fastq format is valid\n"; exit; }	

	# Exteact scaffold sequence 
	#***/ $read_length => read length of fastq format
	#***/ $scaffold => scaffold sequence path
	my %scaff_seq;	# data structure
			# $scaff_seq{$name}[0] = [$tag, $sequence_upstream];	($tag==0: $sequence_upstream == read_length; $tag > 0: $sequence_upstream < read_length)
			# $scaff_seq{$name}[1] = [$tag, $sequence_downstream];	($tag==0: $sequence_downstream == read_length; $tag > 0: $sequence_downstream < read_length)
			# $scaff_seq{$name}[2] = [$tag, $sequence_upstream_reverse_complement];  ($tag==0: $sequence_upstream_reverse_complement == read_length; $tag > 0: $sequence_upstream_reverse_complement < read_length)
			# $scaff_seq{$name}[3] = [$tag, $sequence_downstream_reverse_complement];  ($tag==0: $sequence_downstream_reverse_complement == read_length; $tag > 0: $sequence_downstream_reverse_complement < read_length)
	Check::extract_scaffold($read_length, $scaffold, \%scaff_seq);

	my %uniq_break; # record whether the value (sequence) in %scaff_seq is unique;
	if ( %scaff_seq ) {
		foreach my $name ( keys %scaff_seq ) { # $name => name of breakpoint sequence
			# Check and Rebuild the scafflod sequence
			#***/ $scaff_seq{$name} => see data structure above
			#***/ $tran_seq => cdna sequence data from Ensembl
			#***/ $gene => gene name (symbol) data from Ensembl
			#***/ $geneA and $geneB => GeneA, GeneB given by user
			#***/ $ref_sequence (point to \%scaff_final in module Check::judge_scaffold)
			#	$ref_sequence->{"GeneA"} = [sequence_of_longest_transcript_of_GeneA, ensembl_id_of_longest_transcript_of_GeneA, breakpoint_sequence_of_GeneA] 
			#	$ref_sequence->{"GeneB"} = [sequence_of_longest_transcript_of_GeneB, ensembl_id_of_longest_transcript_of_GeneB, breakpoint_sequence_of_GeneB] 

			print "\n###########################################################################################\n";
			print "# Step 1: Check and build GeneA-GeneB-scaffold sequence for realigning of breakpoint $name #\n";
			print "############################################################################################\n";

			my $tran_seq; # set annotation input path
			my $tag_scaff_geneA; my $tag_scaff_geneB; my $ref_sequence;
			if ( $trans_ref eq "ensembl" ) {
				$tran_seq = "$input/ensembl_transcript.fa"; # cdna mRNA sequences from ensmbl resource
				($tag_scaff_geneA, $tag_scaff_geneB, $ref_sequence) = Check::judge_scaffold($scaff_seq{$name}, $tran_seq, $gene, $geneA, $geneB, "ensembl", $user_ref);
			} elsif ( $trans_ref eq "gencode" ) {
				$tran_seq = "$input/gencode_transcript.fa"; # cdna mRNA sequences from Gencode resource
				($tag_scaff_geneA, $tag_scaff_geneB, $ref_sequence) = Check::judge_scaffold($scaff_seq{$name}, $tran_seq, $gene, $geneA, $geneB, "gencode", $user_ref);
			} elsif ( $trans_ref eq "ucsc" ) {
				$tran_seq = "$input/ucsc_transcript.fa"; # cdna mRNA sequences from ucsc resource
				my $ucsc_name = "$input/ucsc_refGene.txt"; # UCSC gene name and transcript id annotation
				($tag_scaff_geneA, $tag_scaff_geneB, $ref_sequence) = Check::judge_scaffold($scaff_seq{$name}, $tran_seq, $gene, $geneA, $geneB, $ucsc_name, $user_ref);
			}
			my $breakpoint_scaffold; # $breakpoint_scaffold => breakpoint position in scaffold sequence (xxxxx|yyyyy: the position of the first y)

			#assmble the cDNA and scaffold sequence
			if ( $tag_scaff_geneA == 1 ) {
				if ( $tag_scaff_geneB == 1 ) {
					# judge whether the breakpoint sequence is duplicate or unique
					my $sequence = $ref_sequence->{$geneA}[2].$ref_sequence->{$geneB}[2];
					if ( exists($uniq_break{$sequence}) ) { print "Step 1: $name has a duplicate breakpoint sequence\n"; next; } else { $uniq_break{$sequence} = 1; } # if it's duplicate one, go back the loop

					# create a folder named as the given breakpoint (if breakpoint sequence not duplicate one)
					if ( -e "$output/$name" ) {
						`rm -r $output/$name/*`; print "\nWRANINGS: $output/$name exists, and the old file will be re-written\n";
						`mkdir $output/$name/tmp`;
					} else {
						`mkdir $output/$name`;
						`mkdir $output/$name/tmp`;
					}

					# the start position of breakpoint site in the sequence
					$breakpoint_scaffold = length($ref_sequence->{$geneA}[2])+1;
					# output the GeneA, GeneB and breakpoint together as a fasta file (reference for running alignment)
					my $transA = $ref_sequence->{$geneA}[1]; my $transB = $ref_sequence->{$geneB}[1]; # use transcript Ensembl_id (the longest one) as file name
					open (OUT, ">$output/$name/scaffold_${transA}_${breakpoint_scaffold}_${transB}_${read_length}_seq.fa") || die "Step 1: Cannot open this path for scaffold sequence:$!\n";
					print OUT ">$geneA\n";
					print OUT "$ref_sequence->{$geneA}[0]\n";
					print OUT ">$geneB\n";
					print OUT "$ref_sequence->{$geneB}[0]\n";
					print OUT ">scaffold\n"; # Pls note: always use "scaffold" as the name
					print OUT "$ref_sequence->{$geneA}[2]", "$ref_sequence->{$geneB}[2]\n";
					close OUT;
					print "Step 1: The position of breakpoint split in the sequence is at $breakpoint_scaffold, and continues the next step\n";
				} else {
					print "Step 1: Downstream of scaffold does not match the transcript of $geneB and stop here for $name, or use user-defined sequence\n"; next;
				}
			} else {
				if ( $tag_scaff_geneB == 1 ) {
					print "Step 1: Upstream of scaffold does not match the transcript of $geneA and stop here for $name, or use user-defined sequence\n"; next;
				} else {
					print "Step 1: Upstream and downstream of scaffold do not match the transcripts of $geneA and $geneB and stop here for $name, or use user-defined sequence\n"; next;
				}
			}

			print "###############################################################################\n";
			print "# Step 2-1: Align reads to GeneA-GeneB-scaffold sequence for breakpoint $name #\n";
			print "# (no-splicing aligning strategy applied)                                     #\n";
			print "###############################################################################\n";
			#1.1 build index of GeneA-GeneB-scaffold
			my $build_flag = `hisat2-build -p $cpus -f $output/$name/scaffold_*_seq.fa $output/$name/tmp/index`; # name of index file
			if ( $build_flag ) { # whether build success
			} else {
				print "Step 2-1: Build hisat2 index for GeneA-GeneB-scaffold($name) fails\n"; exit;
			}
	
			#1.2 run mapping to GeneA-GeneB-scaffold
			my $hisat_flag = system("hisat2 -p $cpus --no-unal --no-spliced-alignment --no-softclip -x $output/$name/tmp/index -q -1 $fastq_1 -2 $fastq_2 -S $output/$name/tmp/hisats_noclip.sam"); # name of sam file
			$hisat_flag == 0 or die "Step 2-1 Hisat2 mapping to GeneA-GeneB-scaffold fails\n";
			#1.3 transform sam format to sam format using samtools
			if ( -e "$output/$name/tmp/hisats_noclip.sam" ) {
				my $sam_flag = system("samtools view $output/$name/tmp/hisats_noclip.sam -bS > $output/$name/tmp/hisats_noclip.bam");
				$sam_flag == 0 or die "samtools does not work, please check the setting of samtools\n";
				`samtools sort $output/$name/tmp/hisats_noclip.bam -o $output/$name/hisats_noclip.sorted.bam`;
				`samtools index $output/$name/hisats_noclip.sorted.bam`;
				`rm $output/$name/tmp/hisats_noclip.bam`;
			} else {
				print "Step 2-1: Hisat2 map to GeneA-GeneB-scaffold($name): no reads mapping\n"; exit;
			}
		
			print "##################################################################################################\n";
			print "# Step 2-2: Target discordant/singlton discordant reads from the sam file (for breakpoint $name) #\n";
			print "##################################################################################################\n";
			my %spanning = (); #data structure for spanning read
			Extract::Spanning(\%spanning, $geneA, $geneB, "$output/$name/");

			my %discordant = (); #data structure for discordant split read;
			Extract::Discordant(\%discordant, $anchor, $breakpoint_scaffold, $read_length, "$output/$name/", $geneA, $geneB, \%spanning);

			my %singlton = (); #data structure for singlton split read;
			Extract::Singlton(\%singlton, $anchor, $breakpoint_scaffold, $read_length, "$output/$name/");

			print "#########################################################################################\n";
			print "# Step 2-3: Extract discordant and singlton split read from raw fastq file              #\n";
			print "# (parse fastq file, grep the reads that could be aligned to breakpoint sequence $name) #\n";
			print "#########################################################################################\n";
			if ( %discordant ) {
				First_output::grep_discordant(\%discordant, $fastq_1, $fastq_2, "$output/$name/", $header_sep);
			} else {
				print "Step 2-3: No discordant split reads mapped the scaffold\n";
			}

			if ( %singlton ) {
				First_output::grep_singlton(\%singlton, $fastq_1, $fastq_2, "$output/$name/", $header_sep);
			} else {
				print "Step 2-3: no singlton split reads mapped the scaffold\n";
			}

			if ( %spanning ) {
				First_output::grep_spanning(\%spanning, $fastq_1, $fastq_2, "$output/$name/", $header_sep);
			} else {
				print "Step 2-3: no spanning reads mapped the scaffold\n";
			}

			if (! -e "$output/$name/read_mapped_info" ) { `cat <> $output/$name/read_mapped_info`; } # if read_mapped_info not presence, create one with empty
			if (! -e "$output/$name/discordant_split_1.txt" ) { `cat <> $output/$name/discordant_split_1.txt`; `cat <> $output/$name/discordant_split_2.txt`; } # if discordant_split_* not presence, create ones with empty
			if (! -e "$output/$name/singlton_split_1.txt" ) { `cat <> $output/$name/singlton_split_1.txt`; `cat <> $output/$name/singlton_split_2.txt`; } # if singlton_split_* not presence, create ones with empty
			if (! -e "$output/$name/spanning_1.txt" ) { `cat <> $output/$name/spanning_1.txt`; `cat <> $output/$name/spanning_2.txt`; } # if spanning_* not presence, create ones with empty

			print "####################################################################################\n";
			print "# step 3: Realign the select reads to the whole genome (test for specificity)      #\n";
			print "# (Splicing aligning model applied) for breakpoint $name aligning 		  #\n";
			print "####################################################################################\n";
			#Filtering using genome realignment strategy
			if (! -z "$output/$name/read_mapped_info" ) { #whether read_mapped_info is zero
				my %read; #reload the discordant/singlton split reads info file in the step3
					#*****/ %read data structure /
					#$read{$type}{$name}: type = "discordant_split | singlton_split | spanning"; name = read name in fastq file (e.g. "UNC12-SN629:390:C5AAMACXX:2:1306:18427:11314")
					#$read{$type}{$name}[0]: the informaton of first-end read 
					#$read{$type}{$name}[0][0]: mapped GeneA|GeneB|scaffold status (e.g. "scaffold" or "scaffold(reverse_complement)")
					#$read{$type}{$name}[0][1]: orignal fastq sequence of first-end
					#$read{$type}{$name}[0][2]: reverse complementary fastq sequence of first-end
					#$read{$type}{$name}[0][3]: quality info of first-end
					#$read{$type}{$name}[1]: the informaton of second-end read 
					#$read{$type}{$name}[1][0]: mapped GeneA|GeneB|scaffold status (e.g. "scaffold" or "scaffold(reverse_complement)")
					#$read{$type}{$name}[1][1]: orignal fastq sequence of second-end
					#$read{$type}{$name}[1][2]: reverse complementary fastq sequence of second-end
					#$read{$type}{$name}[1][3]: quality info of second-end

        			open (IN, "cut -s -f1-4 $output/$name/read_mapped_info |") || die "Step 3-0: wrong read output for read_mapped_info:$!\n";
        			while ( <IN> ) {
                			chomp $_; next if ( $_ =~/wrong/ );
					my ($type, $name, $first_hit, $second_hit) = (split /\t/, $_)[0,1,2,3];
                			$read{$type}{$name}[0][0] = $first_hit; # name
                			$read{$type}{$name}[1][0] = $second_hit; # name
        			}
        			close IN;

				##################################
				# print the final spanning reads #
				##################################
				#load the gene structure	
				my @partner; # data structure -- collect cDNA sequences of transcriptome
				my ($geneA_pos_tag, $geneB_pos_tag) = Fusion_gene::cDNA(\@partner, $gene, $geneA, $geneB); 

				my %multiple = (); # collect multiple hit reads for discordant_split / singlton_split / spanning
				my $type = 'spanning';
				if ( exists($read{$type}) ) {
					my $spa_flag = system("hisat2 -p $cpus --no-unal --secondary -k 600 -x $genome_index -q -U $output/$name/${type}_1.txt,$output/$name/${type}_2.txt -S $output/$name/tmp/${type}_sec.sam");
					$spa_flag == 0 or die "Step 3-1 Hisat2 maps spanning_1 and spanning_2 to genome fails\n";
					Genome_align::discordant_specif($read{$type}, $multiple{$type}, "$output/$name/", \@partner, "spanning", $header_sep, $geneA_pos_tag, $geneB_pos_tag);
					if ( %{$read{$type}} ) {
						open (OUT1, ">$output/$name/final_spanning_1.txt") || die "Step 3-1: cannot output the final 1st spanning:$!\n";
						open (OUT2, ">$output/$name/final_spanning_2.txt") || die "Step 3-1: cannot output the final 2nd spanning:$!\n";
						open (OUT, ">$output/$name/final_read_mapped_info") || die "Step 3-1: cannot output the final spanning read mapped summary:$!\n";
						foreach my $name ( keys %{$read{$type}} ) {
							print OUT1 "$read{$type}{$name}[0][3]\n";
							print OUT2 "$read{$type}{$name}[1][3]\n";
							print OUT "$type\t$name\t$read{$type}{$name}[0][0]\t$read{$type}{$name}[1][0]\n";
						}
						close OUT1;
						close OUT2;
						close OUT;
					} else {
						`cat <> $output/$name/final_read_mapped_info`; `cat <> $output/$name/final_spanning_1.txt`; `cat <> $output/$name/final_spanning_2.txt`;
					}
				} else { # spanning_* not present; create final_spanning_* as zero
					`cat <> $output/$name/final_read_mapped_info`; `cat <> $output/$name/final_spanning_1.txt`; `cat <> $output/$name/final_spanning_2.txt`;
				}
		
				###############################
				# print the final split reads #
				###############################
				$type = 'discordant_split';
				if ( exists($read{$type}) ) {
					my $dis_flag = system("hisat2 -p $cpus --no-unal --no-softclip --secondary -k 600 -x $genome_index -q -U $output/$name/${type}_1.txt,$output/$name/${type}_2.txt -S $output/$name/tmp/${type}_sec.sam");
					$dis_flag == 0 or die "Step 3-2 Hisat2 maps discordant_split_1 and discordant_split_2 to genome fails\n";
					Genome_align::discordant_specif($read{$type}, $multiple{$type}, "$output/$name/", \@partner, "discordant_split", $header_sep, $geneA_pos_tag, $geneB_pos_tag);
					if ( %{$read{$type}} ) {
						open (OUT1, ">$output/$name/final_split_1.txt") || die "Step 3-2: cannot output the final 1st discordant split:$!\n";
						open (OUT2, ">$output/$name/final_split_2.txt") || die "Step 3-2: cannot output the final 2nd discordant split:$!\n";
						open (OUT, ">>$output/$name/final_read_mapped_info") || die "Step 3-2: cannot output the final discordant split read mapped summary:$!\n";
						foreach my $name ( keys %{$read{$type}} ) {
							print OUT1 "$read{$type}{$name}[0][3]\n";
							print OUT2 "$read{$type}{$name}[1][3]\n";
							print OUT "$type\t$name\t$read{$type}{$name}[0][0]\t$read{$type}{$name}[1][0]\n";
						}
						close OUT1;
						close OUT2;
						close OUT;
					} else {
						`cat <> $output/$name/final_split_1.txt`; `cat <> $output/$name/final_split_2.txt`;
					}
				} else {
					`cat <> $output/$name/final_split_1.txt`; `cat <> $output/$name/final_split_2.txt`;
				}
	
				$type = 'singlton_split';
				if ( exists($read{$type}) ) {
					my $sing_flag = system("hisat2 -p $cpus --no-unal --no-softclip --secondary -k 600 -x $genome_index -q -U $output/$name/${type}_1.txt,$output/$name/${type}_2.txt -S $output/$name/tmp/${type}_sec.sam");
					$sing_flag == 0 or die "Step 3-3 Hisat2 maps singlton_split_1 and singlton_split_2 to genome fails\n";
					Genome_align::singlton_specif($read{$type}, $multiple{$type}, "$output/$name", \@partner, "singlton_split", $header_sep, $geneA_pos_tag, $geneB_pos_tag);
					if ( %{$read{$type}} ) {
						open (OUT1, ">>$output/$name/final_split_1.txt") || die "Step 3-3: cannot output the final 1st singlton split:$!\n";
						open (OUT2, ">>$output/$name/final_split_2.txt") || die "Step 3-3: cannot output the final 2nd singlton split:$!\n";
						open (OUT, ">>$output/$name/final_read_mapped_info") || die "Step 3-3: cannot output the final singlton split read mapped summary:$!\n";
						foreach my $name ( keys %{$read{$type}} ) {
							print OUT1 "$read{$type}{$name}[0][3]\n";
							print OUT2 "$read{$type}{$name}[1][3]\n";
							if ( $read{$type}{$name}[0][0] eq "NULL" or $read{$type}{$name}[1][0] eq "NULL" ) {
								print OUT "$type\t$name\t$read{$type}{$name}[0][0]\t$read{$type}{$name}[1][0]\n";
							} else { # singlton to discordant split read
								print OUT "singlton_to_discordant_split\t$name\t$read{$type}{$name}[0][0]\t$read{$type}{$name}[1][0]\n";
							} 
						}
						close OUT1;
						close OUT2;
						close OUT;
					}
				}

				# summarize multiple / unique mapping read
				my $Dis_A = `grep 'discordant_split' $output/$name/read_mapped_info | cut -s -f1,3,4 | grep -c $geneA`; chomp $Dis_A;
				if ( $Dis_A eq "" ) { $Dis_A = 0; }
				my $Dis_B = `grep 'discordant_split' $output/$name/read_mapped_info | cut -s -f1,3,4 | grep -c $geneB`; chomp $Dis_B;
				if ( $Dis_B eq "" ) { $Dis_B = 0; }
				my $Spa_AB = `grep -c -P 'spanning\t' $output/$name/read_mapped_info`; chomp $Spa_AB;
				if ( $Spa_AB eq "" ) { $Spa_AB = 0; }

				my $Dis_A_F = `grep 'discordant_split' $output/$name/final_read_mapped_info | cut -s -f1,3,4 | grep -c $geneA`; chomp $Dis_A_F;
				if ( $Dis_A_F eq "" ) { $Dis_A_F = 0; }
				my $Dis_B_F = `grep 'discordant_split' $output/$name/final_read_mapped_info | cut -s -f1,3,4 | grep -c $geneB`; chomp $Dis_B_F;
				if ( $Dis_B_F eq "" ) { $Dis_B_F = 0; }
				my $Spa_AB_F = `grep -c -P 'spanning\t' $output/$name/final_read_mapped_info`; chomp $Spa_AB_F;
				if ( $Spa_AB_F eq "" ) { $Spa_AB_F = 0; }

				my $Per_Dis_A; my $Per_Dis_B; my $Per_Spa_AB;
				if ( $Dis_A == 0 ) { $Per_Dis_A = 0; } else { $Per_Dis_A = $Dis_A_F/$Dis_A; }
				if ( $Dis_B == 0 ) { $Per_Dis_B = 0; } else { $Per_Dis_B = $Dis_B_F/$Dis_B; }
				if ( $Spa_AB == 0 ) { $Per_Spa_AB = 0; } else { $Per_Spa_AB = $Spa_AB_F/$Spa_AB; }

				open (SPT, ">$output/$name/summary_read_mapping_support.txt") || die "Step 3-4: cannot output the summary of read mapping support:$!\n";
				print SPT "# summary of read mapping support #\n";
				print SPT "===================================\n";
				print SPT "* Num (proportion) of discordant split reads unique mapping to $geneA and scaffold: $Dis_A_F / $Dis_A (", $Per_Dis_A, ")\n";
				print SPT "* Num (proportion) of discordant split reads unique mapping to $geneB and scaffold: $Dis_B_F / $Dis_B (", $Per_Dis_B, ")\n";
				print SPT "* Num (proportion) of spanning read pairs unique mapping to $geneA and $geneB: $Spa_AB_F / $Spa_AB (", $Per_Spa_AB, ")\n\n";
				close SPT;

				if (! -z "$output/$name/final_read_mapped_info" ) { # final_read_mapped_info is not zero => put spanning and split reads together for running alignment
					if (! -z "$output/$name/final_spanning_1.txt" ) { # final_spanning_read is not zero => running final spanning read alignment
						my $final_flag = system("hisat2 -p $cpus --no-unal --no-spliced-alignment --no-softclip -x $output/$name/tmp/index -q -1 $output/$name/final_spanning_1.txt -2 $output/$name/final_spanning_2.txt -S $output/$name/tmp/final_spanning_noclip.sam"); # name of sam file (spanning)
        					$final_flag == 0 or die "Step 3-4 Hisat2 map to GeneA-GeneB-scaffold fails (spanning)\n";
						# 4.4 transform sam format to sam format using samtools (spanning)
						if ( -e "$output/$name/tmp/final_spanning_noclip.sam" ) {
							`samtools view $output/$name/tmp/final_spanning_noclip.sam -bS > $output/$name/tmp/final_spanning_noclip.bam`;
							`samtools sort $output/$name/tmp/final_spanning_noclip.bam -o $output/$name/final_spanning_noclip.sorted.bam`;
							`samtools index $output/$name/final_spanning_noclip.sorted.bam`;
							`rm $output/$name/tmp/final_spanning_noclip.bam`;
							`rm $output/$name/tmp/final_spanning_noclip.sam`;
						}
					}
					if ( ! -z "$output/$name/final_split_1.txt" ) { # final_split_read is not zero => running final split read alignment
						my $final_flag = system("hisat2 -p $cpus --no-unal --no-spliced-alignment --no-softclip -x $output/$name/tmp/index -q -1 $output/$name/final_split_1.txt -2 $output/$name/final_split_2.txt -S $output/$name/tmp/final_split_noclip.sam"); # name of sam file (split)
						$final_flag == 0 or die "Step 3-4 Hisat2 map to GeneA-GeneB-scaffold fails (split)\n";
						# 4.4 transform sam format to sam format using samtools (split)
						if ( -e "$output/$name/tmp/final_split_noclip.sam" ) {
							`samtools view $output/$name/tmp/final_split_noclip.sam -bS > $output/$name/tmp/final_split_noclip.bam`;
							`samtools sort $output/$name/tmp/final_split_noclip.bam -o $output/$name/final_split_noclip.sorted.bam`;
							`samtools index $output/$name/final_split_noclip.sorted.bam`;
							`rm $output/$name/tmp/final_split_noclip.bam`;
							`rm $output/$name/tmp/final_split_noclip.sam`;
						}
					}
				}
			} else {
				`cat <> $output/$name/final_read_mapped_info`; 
				`cat <> $output/$name/final_spanning_1.txt`; `cat <> $output/$name/final_spanning_2.txt`;
				`cat <> $output/$name/final_split_1.txt`; `cat <> $output/$name/final_split_2.txt`;
				`echo '# summary of read mapping support #' > $output/$name/summary_read_mapping_support.txt`;
				`echo '===================================' >> $output/$name/summary_read_mapping_support.txt`;
				`echo "* Num (proportion) of discordant split reads unique mapping to $geneA and scaffold: 0 / 0 (0)" >> $output/$name/summary_read_mapping_support.txt`;
				`echo "* Num (proportion) of discordant split reads unique mapping to $geneB and scaffold: 0 / 0 (0)" >> $output/$name/summary_read_mapping_support.txt`;
				`echo "* Num (proportion) of spanning read pairs unique mapping to $geneA and $geneB: 0 / 0 (0)" >> $output/$name/summary_read_mapping_support.txt`;
				print "Step 3: No spanning | discordant/singlton split reads show successful mapping for breakpoint $name\n";
			}
		}
	}

