#!/usr/bin/perl -w

########################################################################################################################################################
# Evaluate.pl is perl script that concatenates split read across all samples and evaluate the read coverage and distribution that cross the breakpoint 
# Before you start evaluate.pl, you have to run Select_read.pl first
# The script was created in March 10, 2017, and modified in Sept 13, 2017
########################################################################################################################################################
use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Spec;

	# setup option parameters in the script command line
	my @usage;
	push @usage, "Usage: ".basename($0)." [options]\n";
	push @usage, "Evaluate split reads that support a given fusion breakpoint.\n";
	push @usage, "  --help          Displays this information\n";
	push @usage, "  --input         Set input path for running the script: this path is the output path of running Select_read.pl (it should contain all the samples that you are going to evaluate)\n";
	push @usage, "  --output        Set output path for running the script\n";

	my $help;
	my $input;
	my $output;
	GetOptions
	(
		'help'        => \$help,
		'input=s'     => \$input,
		'output=s'    => \$output
	);
	not defined $help or die @usage;
	defined $input or die @usage;
	defined $output or die @usage;

	# setup output folder status
	if ( -e $output ) {
		`rm -r $output/*`; print "\nwarning: $output exists, and the old files will be re-written\n";
	} else {
		`mkdir $output`;
	}

	my %break_num;	# data structure: $break_num{$break}{$seq_length} = [$break_start, $read_length];
			# $break = unique breakpoint name (e.g. alt_1)
			# $seq_length = the length of scaffold sequence
			# [$break_start: the start position of breakpoint in the scaffold sequence; $read_length: the length of pairend reads; $scaffold_name: the scaffold name of one given sample]
	my %summary; # summarize spanning, discordant_split, singlton_split reads across all samples
	if ( -e "$input" ) {	
		open (INPUT, "ls $input |") || die "main $input path is wrong\n";
		while ( <INPUT> ) {
			chomp $_; my $name = $_; # sample name
			if ( -d "$input/$name" ) { # judge whether sample is a folder
				open (BK, "ls $input/$name |") || die "main $input/$name path is wrong\n";
				while ( <BK> ) {
					chomp $_; my $break = $_; print "### $name && $break ###\n"; # name of breakpoint
					# if breakpoint folder not present in output path, please create it.
					if (! -e "$output/$break") { `mkdir $output/$break`; } 
					
					if ( -d "$input/$name/$break" ) { # judge whether breakpoint is folder for input path

						if ( -e "$input/$name/$break/final_split_1.txt" ) { # judge whether final_split_1/_2 file present
							if (! -z "$input/$name/$break/final_split_1.txt" ) { # judge whether final_split_1/_2 is not zero
								`cat $input/$name/$break/final_split_1.txt >> $output/$break/All_sample_split_1.txt`;
								`cat $input/$name/$break/final_split_2.txt >> $output/$break/All_sample_split_2.txt`;
							}
							# manage the breakpoint sequence: extract three important paramters (breakpoint_start, scaffold_length, read_length)
							my $tmp = `ls $input/$name/$break/scaffold_ENST*_ENST*_seq.fa`; chomp $tmp;;
							my $file_name = (split /\//, $tmp)[-1]; # the last element
							my ($break_start, $read_length) = (split /_/, $file_name)[2, 4];
							if ( defined($break_start) && defined($read_length) ) {
								my $sequence = `grep -A1 'scaffold' $tmp`; chomp $sequence; $sequence = (split /\n/, $sequence)[1]; 
								my $seq_length = length($sequence);
								$break_num{$break}{$seq_length} = [$break_start, $read_length, $tmp, $sequence];
								#print "Breakpoint_sequence: $break\t$seq_length\t$break_start, $read_length, $tmp, $sequence\n";
							} else {
								print "scaffold name have been changed, please apply original one: scaffold_ENST*_XX_ENST*_YY_seq.fa\n";
							}
						}

						if ( -e "$input/$name/$break/final_read_mapped_info" ) { # judge whether final_read_mapped_info present
							my $num_span = 0; # number of spanning reads (one mapped to geneA, the other mapped to geneB)
							my $num_split_dis = 0; # number of discordant_split read (one mapped to scaffold, the other mapped to geneA / geneB)
							my $num_split_sig = 0; # number of singlton_split read ( one mapped to scaffold, the other shows no mapping)
							my $num_split_dis_sig = 0; # number of singlton_split => discordant_split reads during genome mapping process

							if (! -z "$input/$name/$break/final_read_mapped_info" ) { # judge whether final_read_mapped_info is not zero
								$num_span = `grep -c -P '^spanning' $input/$name/$break/final_read_mapped_info`; chomp $num_span;
								$num_split_dis = `grep -c -P '^discordant_split' $input/$name/$break/final_read_mapped_info`; chomp $num_split_dis;
								$num_split_dis_sig = `grep -c -P '^singlton_to_discordant_split' $input/$name/$break/final_read_mapped_info`; chomp $num_split_dis_sig;
								$num_split_sig = `grep -c -P '^singlton_split' $input/$name/$break/final_read_mapped_info`; chomp $num_split_sig;
							}
							$summary{$name}{$break} = [$num_split_dis, $num_split_dis_sig, $num_split_sig, $num_span];
							#print "Sample_Break: $name\t$break\t$num_split_dis, $num_split_dis_sig, $num_split_sig, $num_span\n";
						}
					}

				}
				close BK;
			}
		}
		close INPUT;
	}

	# run statistical test for balance mapping of breakpoint scaffold
	my %summary_pvalue; # collect [num_read_upstream, num_read_downstream, pvalue] per one breakpoint scaffold
	foreach my $break ( keys %break_num ) { # $break represent name of unique breakpoint
		my @array = sort {$a <=> $b} keys %{$break_num{$break}}; # Get all the key %{$break_num{$break}} (scaffold sequence length);

		my $upstream = 0; my $downstream = 0; my $pvalue = "NULL";
		if ( scalar(@array) == 1 ) { # only one unique scaffold sequence, indicating samples with the same read length
			# cp the scaffold sequence to the output folder
			`cp $break_num{$break}{$array[0]}[2] $output/$break`;
			if (! -e "$output/$break/All_sample_split_1.txt" ) { `cat <> $output/$break/All_sample_split_1.txt`; `cat <> $output/$break/All_sample_split_1.txt`; } # NULL

			if (! -z "$output/$break/All_sample_split_1.txt" ) {
				my $build_flag = `hisat2-build -p 1 -f $output/$break/scaffold_ENST*_ENST*_seq.fa $output/$break/index`; # build index file
				my $final_flag = system("hisat2 -p 1 --no-unal --no-spliced-alignment --no-softclip -x $output/$break/index -q -1 $output/$break/All_sample_split_1.txt -2 $output/$break/All_sample_split_2.txt -S $output/$break/All_sample_noclip.sam");
				$final_flag == 0 or die "Hisat2 map to GeneA-GeneB-scaffold for all concatenated samples fails\n";

				# transform sam format to sam format using samtools
				`samtools view $output/$break/All_sample_noclip.sam -bS > $output/$break/All_sample_noclip.bam`;
				`samtools sort $output/$break/All_sample_noclip.bam -o $output/$break/All_sample_noclip.sorted.bam`;
				`samtools index $output/$break/All_sample_noclip.sorted.bam`;
				`rm $output/$break/All_sample_noclip.bam`;
				`rm $output/$break/index*`;
				
				# set the three parameters
				# transfer the three parameters: scaffold_length, the start position of breakpoint site, read_length, output_path
				($upstream, $downstream, $pvalue) = &access_distribution($array[0], $break_num{$break}{$array[0]}[0], $break_num{$break}{$array[0]}[1], "$output/$break"); 
			}
		} else { # there are more than one scaffold sequence, indicating samples with various read length
			`cp $break_num{$break}{$array[-1]}[2] $output/$break`; # put the longest scaffold as reference candidate
			if (! -e "$output/$break/All_sample_split_1.txt" ) { `cat <> $output/$break/All_sample_split_1.txt`; `cat <> $output/$break/All_sample_split_1.txt`; } # NULL

			if (! -z "$output/$break/All_sample_split_1.txt" ) {
				my $build_flag = `hisat2-build -p 1 -f $output/$break/scaffold_ENST*_ENST*_seq.fa $output/$break/index`; # build index file
				my $final_flag = system("hisat2 -p 1 --no-unal --no-spliced-alignment --no-softclip -x $output/$break/index -q -1 $output/$break/All_sample_split_1.txt -2 $output/$break/All_sample_split_2.txt -S $output/$break/All_sample_noclip.sam");
				$final_flag == 0 or die "Hisat2 map to GeneA-GeneB-scaffold for all concatenated samples fails\n";

				# transform sam format to sam format using samtools
				`samtools view $output/$break/All_sample_noclip.sam -bS > $output/$break/All_sample_noclip.bam`;
				`samtools sort $output/$break/All_sample_noclip.bam -o $output/$break/All_sample_noclip.sorted.bam`;
				`samtools index $output/$break/All_sample_noclip.sorted.bam`;
				`rm $output/$break/All_sample_noclip.bam`;
				`rm $output/$break/index*`;
				
				# set the three parameters
				# transfer the three parameters: scaffold_length, the start position of breakpoint site, read_length, output_path
				($upstream, $downstream, $pvalue) = &access_distribution($array[-1], $break_num{$break}{$array[-1]}[0], $break_num{$break}{$array[0]}[1], "$output/$break"); 
			}
		}
		open (OUT, ">$output/$break/p.value") || die "cannot output pvalue:$!\n";
		print OUT "Num_mapped_upstream\tNum_mapped_downstream\tp.value(Fisher_test)\n";
		print OUT "$upstream\t$downstream\t$pvalue\n";
		close OUT;
		$summary_pvalue{$break} = [$upstream, $downstream, $pvalue];
	}

	# final print the output #
	open (OUT, ">$output/summary.txt") || die "cannot open output file:$!\n";
	print OUT "Sample";
	foreach my $break ( sort {$a cmp $b} keys %summary_pvalue ) { print OUT "\t$break"; }
	print OUT "\n";
	foreach my $name ( sort {$a cmp $b} keys %summary ) {
		print OUT "$name";
		foreach my $break ( sort {$a cmp $b} keys %{$summary{$name}} ) {
			print OUT "\t$summary{$name}{$break}[0]", "+", "$summary{$name}{$break}[1],", "$summary{$name}{$break}[2];", "$summary{$name}{$break}[3]";
		}
		print OUT "\n";
	}
	print OUT "Statistics_summary(split_reads)";
	foreach my $break ( sort {$a cmp $b} keys %summary_pvalue ) {
		print OUT "\t$summary_pvalue{$break}[0]|$summary_pvalue{$break}[1]: $summary_pvalue{$break}[2]";
	}
	print OUT "\n";

########sub function##########
	sub access_distribution {
		my ($scaffold_length, $break_start, $read_length, $path) = @_;
		my %collect; # collect the read names 
		# Load sam file
		open (IN, "awk -F '\t' -v OFS='\t' '(\$8>0)' $path/All_sample_noclip.sam |") || die "cannot run awk scaffold recruit script:$!\n";
		while ( <IN> ) {
			my ($name, $one_match, $one_pos, $other_match, $other_pos, $seq) = (split /\t/, $_)[0, 2, 3, 6, 7, 9];
			$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g;
			if ( $one_match eq "scaffold" ) {
				push @{$collect{$name}}, [$one_pos, $seq];
			}
		}
		close IN;

		my %ref; # statistics on read mapping bias in scaffold
		foreach my $name ( keys %collect ) {
			if ( scalar(@{$collect{$name}}) == 1 ) {
				my $middle = $collect{$name}[0][0] + length($collect{$name}[0][1])/2; # the middle position of reads
				if ( $middle > $break_start ) {
					push @{$ref{"downstream"}}, $name;
				} elsif ( $middle <= $break_start ) {
					push @{$ref{"upstream"}}, $name;
				}
			} elsif ( scalar(@{$collect{$name}}) == 2 ) {
				if ( $collect{$name}[0][0] < $collect{$name}[1][0] ) { # start point: first element < second element
					my $middle = ($collect{$name}[0][0] + $collect{$name}[1][0] + length($collect{$name}[1][1]))/2;
					if ( $middle > $break_start ) {
						push @{$ref{"downstream"}}, $name;
					} elsif ( $middle <= $break_start ) {
						push @{$ref{"upstream"}}, $name;
					}
				} else {
					my $middle = ($collect{$name}[1][0] + $collect{$name}[0][0] + length($collect{$name}[0][1]))/2;
					if ( $middle > $break_start ) {
						push @{$ref{"downstream"}}, $name;
					} elsif ( $middle <= $break_start ) {
						push @{$ref{"upstream"}}, $name;
					}
				}
			}
		}

		#run statistical anaylses
		my $up = 0; my $down = 0; my $p;
		if ( exists($ref{"upstream"}) ) {
			if ( exists($ref{"downstream"}) ) {
				$up = scalar(@{$ref{"upstream"}});
				$down = scalar(@{$ref{"downstream"}});
			} else {
				$up = scalar(@{$ref{"upstream"}});
			}
		} else {
			if ( exists($ref{"downstream"}) ) {
				$down = scalar(@{$ref{"downstream"}});
			}
		}

		if ( $up == 0 && $down == 0 ) { #if no reads mapped to scaffold, please cancel the statistical test
			return($up, $down, "NA"); # please reture NA for pvalue
		} else {
			`echo "a=$up; b=$down; a1=round(($up+$down)*($break_start/$scaffold_length), digits=0); b1=round(($up+$down)*(1-$break_start/$scaffold_length), digits=0); ab=matrix(c(a,b,a1,b1),nrow=2); c=fisher.test(ab); write.table(c\$ p.value, file='$path/p.value', row.names=F, col.names=F); sessionInfo()" | R --vanilla --slave`;
			my $p = `head -n 1 $path/p.value`; chomp $p;
			return($up, $down, $p);
		}
	}

