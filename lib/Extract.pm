package Extract;
use strict;
use warnings;

	sub Spanning {
		my ($dis_ref, $geneA, $geneB, $sam_file, $aligner) = @_;

		# select mapped spanning reads (one mapping to GeneA/GeneB; the other mapping to GeneB/GeneA)
		my %spanning_read; # collect the mapped spanning reads
		open (IN, "awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7!=\"=\" && \$5>=0 && \$8>0)' $sam_file/tmp/noclip.sam | ") || die "Step2-2 spanning read mapping: cannot run awk 1:$!\n";
		while ( <IN> ) {
			chomp $_; #
			my ($name, $flag, $one_match, $one_pos, $value, $other_match, $other_pos, $seq, $quality) = (split /\t/, $_)[0,1,2,3,4,6,7,9,10];
			$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # trimmed header
			my $string = ""; # record mismatching information
			if ( $_ =~/MD\:Z\:([\w\^]+)/ ) { 
				$string = $1;
			} else {
				print "Step2-2 spanning read mapping: $name mismatch string contains unexpected symbols\n";
			}
			$string =~s/\^[A-Za-z]+//g; $string =~s/[0-9]+//g; my $num_mis = length($string); # the number of mismatch nucleatides (SNPs+Indels) in read mapping

			if ( $num_mis < 4 ) { # the default setting -- maximum allowed mismatching number
				push @{$spanning_read{$name}}, [$one_pos, "spanning", $one_match, $seq, $other_match, $flag, $quality];
			}
		}
		close IN;
			
		if ( $aligner eq 'hisat2' ) {
			foreach my $name ( keys %spanning_read ) {
				if ( scalar(@{$spanning_read{$name}}) == 2 ) {
					if ( $spanning_read{$name}[0][2] eq $geneA && $spanning_read{$name}[0][4] eq $geneB ) {
						if ( $spanning_read{$name}[1][2] eq $geneB && $spanning_read{$name}[1][4] eq $geneA ) {
							$dis_ref->{$name} = $spanning_read{$name};
						}
					}
					if ( $spanning_read{$name}[0][2] eq $geneB && $spanning_read{$name}[0][4] eq $geneA ) {
						if ( $spanning_read{$name}[1][2] eq $geneA && $spanning_read{$name}[1][4] eq $geneB ) {
							$dis_ref->{$name} = $spanning_read{$name};
						}
					}
				} elsif ( scalar(@{$spanning_read{$name}}) == 3 ) {
					my @uniq_name = do { my %name_count; grep { !$name_count{$_}++ } ($spanning_read{$name}[0][2], $spanning_read{$name}[1][2], $spanning_read{$name}[2][2]) };
					if ( scalar(@uniq_name) == 3 ) {
						my $flag_count = 0; # select the read with correct flag values 
						foreach my $value ( @{$spanning_read{$name}} ) {
							if ( $value->[5] == 81 or $value->[5] == 161 or $value->[5] == 97 or $value->[5] == 145 or $value->[5] == 65 or $value->[5] == 129 or $value->[5] == 113 or $value->[5] == 177 ) { # acceptable flag id
								$flag_count++;
							}
						}
						if ( $flag_count == 2 ) {
							foreach my $id ( @{$spanning_read{$name}} ) {
								if ( $id->[2] eq $geneA ) {
									push @{$dis_ref->{$name}}, $id;
								} elsif ( $id->[2] eq $geneB ) {
									push @{$dis_ref->{$name}}, $id;
								}
							}
						} else {
							print "Step2-2 spanning read mapping from hisat2: $name shows incorrect spanning reads (flag)\n";
						}
					} else {
						print "Step2-2 spanning read mapping from hisat2: $name shows incorrect spanning reads\n";
					}
				}
			}
		} elsif ( $aligner eq 'star' ) {
			# select multiple mapped read in *Aligned.out.sam file for STAR
			my %multi_read; # specific for star
			open (IN, "awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7==\"*\" && \$5<=1 && \$8==0)' $sam_file/tmp/noclipAligned.out.sam | ") || die "Step2-2 select read mapping for star under no-splicing mode: cannot run awk 1:$!\n";
			while ( <IN> ) {
				chomp $_; #
				my ($name, $flag, $one_match, $one_pos, $value, $other_match, $other_pos, $seq, $quality) = (split /\t/, $_)[0,1,2,3,4,6,7,9,10];
				$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # trimmed header
				push @{$multi_read{$name}}, [$one_pos, "spanning", $one_match, $seq, $other_match, $flag, $quality];
			}
			close IN;

			foreach my $name ( keys %multi_read ) {
				if ( exists($spanning_read{$name}) ) { # if multi_mapping reads present in %spanning_read
					next if ( scalar(@{$spanning_read{$name}}) == 1 );
					if ( scalar(@{$multi_read{$name}}) == 3 ) {
						my %hash_tmp; my @judge_array = grep { ++$hash_tmp{$_} == 1 } ($multi_read{$name}[0][2], $multi_read{$name}[1][2], $multi_read{$name}[2][2]);
						if ( scalar(@judge_array) == 3 ) { # pair-ended reads can map to GeneA, GeneB and scaffold, respectively
							if ( scalar(@{$spanning_read{$name}}) == 2 ) {
								if ( $spanning_read{$name}[0][2] eq $geneA && $spanning_read{$name}[0][4] eq "scaffold" ) {  # first geneA, second scaffold
									$spanning_read{$name}[0][4] = "$geneB"; # change $other_match of first read to geneB
									$spanning_read{$name}[1][2] = $geneB; # change $one_match of second read to that of geneB
									$spanning_read{$name}[1][4] = $geneA; # change $other_match of second read to geneA
									if ( $multi_read{$name}[0][2] eq "$geneB" ) {
										$spanning_read{$name}[1][0] = $multi_read{$name}[0][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[1][2] eq "$geneB" ) {
										$spanning_read{$name}[1][0] = $multi_read{$name}[1][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[2][2] eq "$geneB" ) {
										$spanning_read{$name}[1][0] = $multi_read{$name}[2][0];  # change $one_pos of second read to that of geneB
									}
								} elsif ( $spanning_read{$name}[0][2] eq "scaffold" && $spanning_read{$name}[0][4] eq $geneA ) { # first scaffold, second geneA
									$spanning_read{$name}[1][4] = "$geneB"; # change $other_match of first read to geneB
									$spanning_read{$name}[0][2] = $geneB; # change $one_match of second read to that of geneB
									$spanning_read{$name}[0][4] = $geneA; # change $other_match of second read to geneA
									if ( $multi_read{$name}[0][2] eq "$geneB" ) { 
										$spanning_read{$name}[0][0] = $multi_read{$name}[0][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[1][2] eq "$geneB" ) {
										$spanning_read{$name}[0][0] = $multi_read{$name}[1][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[2][2] eq "$geneB" ) {
										$spanning_read{$name}[0][0] = $multi_read{$name}[2][0];  # change $one_pos of second read to that of geneB
									}
								}

								if ( $spanning_read{$name}[0][2] eq $geneB && $spanning_read{$name}[0][4] eq "scaffold" ) { # first geneB, second scaffold
									$spanning_read{$name}[0][4] = $geneA; # change $other_match of first read to geneA
									$spanning_read{$name}[1][2] = $geneA; # change $one_match of second read to that of geneA
									$spanning_read{$name}[1][4] = $geneB; # change $other_match of second read to geneB
									if ( $multi_read{$name}[0][2] eq "$geneA" ) { 
										$spanning_read{$name}[1][0] = $multi_read{$name}[0][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[1][2] eq "$geneA" ) {
										$spanning_read{$name}[1][0] = $multi_read{$name}[1][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[2][2] eq "$geneA" ) {
										$spanning_read{$name}[1][0] = $multi_read{$name}[2][0];  # change $one_pos of second read to that of geneB
									}
								} elsif ( $spanning_read{$name}[0][2] eq "scaffold" && $spanning_read{$name}[0][4] eq $geneB ) { # first scaffold, second geneB
									$spanning_read{$name}[1][4] = "$geneA"; # change $other_match of first read to geneB
									$spanning_read{$name}[0][2] = $geneA; # change $one_match of second read to that of geneB
									$spanning_read{$name}[0][4] = $geneB; # change $other_match of second read to geneB
									if ( $multi_read{$name}[0][2] eq "$geneA" ) {
										$spanning_read{$name}[0][0] = $multi_read{$name}[0][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[1][2] eq "$geneA" ) {
										$spanning_read{$name}[0][0] = $multi_read{$name}[1][0];  # change $one_pos of second read to that of geneB
									} elsif ( $multi_read{$name}[2][2] eq "$geneA" ) {
										$spanning_read{$name}[0][0] = $multi_read{$name}[2][0];  # change $one_pos of second read to that of geneB
									}
								}
							}
						} else {
							delete($spanning_read{$name});
						}
					} elsif ( scalar(@{$multi_read{$name}}) > 3 ) {
						delete($spanning_read{$name});
					}
				} else { # if multi_mapping reads not present in %spanning_read
					if ( scalar(@{$multi_read{$name}}) == 3 ) { # add this condition to %spanning_read
						my %hash_tmp; my @judge_array = grep { ++$hash_tmp{$_} == 1 } ($multi_read{$name}[0][2], $multi_read{$name}[1][2], $multi_read{$name}[2][2]);
						if ( scalar(@judge_array) == 3 ) { # pair-ended reads can map to GeneA, GeneB and scaffold, respectively
							if ( $multi_read{$name}[0][2] eq "scaffold" ) {
								$spanning_read{$name}[0] = [$multi_read{$name}[1][0], "spanning", $multi_read{$name}[1][2], $multi_read{$name}[1][3], $multi_read{$name}[2][2], 81, $multi_read{$name}[1][6]];
								$spanning_read{$name}[1] = [$multi_read{$name}[2][0], "spanning", $multi_read{$name}[2][2], $multi_read{$name}[2][3], $multi_read{$name}[1][2], 161, $multi_read{$name}[2][6]];
							} elsif ( $multi_read{$name}[1][2] eq "scaffold" ) {
								$spanning_read{$name}[0] = [$multi_read{$name}[0][0], "spanning", $multi_read{$name}[0][2], $multi_read{$name}[0][3], $multi_read{$name}[2][2], 81, $multi_read{$name}[0][6]];
								$spanning_read{$name}[1] = [$multi_read{$name}[2][0], "spanning", $multi_read{$name}[2][2], $multi_read{$name}[2][3], $multi_read{$name}[0][2], 161, $multi_read{$name}[2][6]];
							} elsif ( $multi_read{$name}[2][2] eq "scaffold" ) {
								$spanning_read{$name}[0] = [$multi_read{$name}[0][0], "spanning", $multi_read{$name}[0][2], $multi_read{$name}[0][3], $multi_read{$name}[1][2], 81, $multi_read{$name}[0][6]];
								$spanning_read{$name}[1] = [$multi_read{$name}[1][0], "spanning", $multi_read{$name}[1][2], $multi_read{$name}[1][3], $multi_read{$name}[0][2], 161, $multi_read{$name}[1][6]];
							}
						}
					}
				}
			}

			foreach my $name ( keys %spanning_read ) {	
				if ( scalar(@{$spanning_read{$name}}) == 2 ) {
					if ( ($spanning_read{$name}[0][2] eq $geneA && $spanning_read{$name}[0][4] eq $geneB) || ($spanning_read{$name}[0][2] eq $geneB && $spanning_read{$name}[0][4] eq $geneA) ) {
						$dis_ref->{$name} = $spanning_read{$name};
					}
				} else {
					print "Step2-2 spanning read mapping for star: $name shows incorrect spanning reads, multiple mapped sites\n";	
				}
			}
		}	
	}
		
	sub Discordant {
		my ($dis_ref, $dis_to_sing_ref, $archor, $breakpoint_scaffold, $read_length, $sam_file, $geneA, $geneB, $dis_ref_span, $aligner) = @_;

		my %discordant_split; # select mapped discordant split reads (one mapping to GeneA/GeneB; the other mapping to scaffold)
		my %multi_read; # select multiple mapped read in *Aligned.out.sam file for STAR
		if ( $aligner eq "hisat2" ) {
			#*****/ 1. mapping quality >=40; ($3!=$7; $7!="=" -- discordant: one mapped to scaffold; the other mapped to GeneA/GeneB) for hisat2/
			open (IN, "awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7!=\"=\" && \$5>=40 && \$8>0)' $sam_file/tmp/noclip.sam | grep 'caffold' |") || die "Step2-2 discordant-split read mapping for hisat2: cannot run awk 1:$!\n";
		} elsif ( $aligner eq "star" ) {
			#*****/ 1. mapping quality <=1; ($3!=$7; $7=="*" -- multiple mapping: at least one of paired has unspecific loci mapped) for star /
			open (IN, "awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7==\"*\" && \$5<=1 && \$8==0)' $sam_file/tmp/noclipAligned.out.sam | ") || die "Step2-2 select read mapping for star under no-splicing mode: cannot run awk 1:$!\n";
			while ( <IN> ) {
				chomp $_; #
				my ($name, $flag, $one_match, $one_pos, $value, $other_match, $other_pos, $seq, $quality) = (split /\t/, $_)[0,1,2,3,4,6,7,9,10];
				$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # trimme header
				push @{$multi_read{$name}}, [$one_pos, "spanning", $one_match, $seq, $other_match, $flag, $quality];
			}
			close IN;
			#*****/ 1. mapping quality >=3; ($3!=$7; $7!="=" -- discordant: one mapped to scaffold; the other mapped to GeneA/GeneB) for star /
			open (IN, "awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7!=\"=\" && \$5>=3 && \$8>0)' $sam_file/tmp/noclip.sam | grep 'caffold' |") || die "Step2-2 discordant-split read mapping for star: cannot run awk 1:$!\n";
		}

		while ( <IN> ) {
			chomp $_; #
			my ($name, $flag, $one_match, $one_pos, $other_match, $other_pos, $seq, $quality) = (split /\t/, $_)[0,1,2,3,6,7,9,10];
			$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # trimmed header
			my $string = ""; # record mismatching information
			if ( $_ =~/MD\:Z\:([\w\^]+)/ ) { 
				$string = $1;
			} else {
				print "Step2-2 discordant-split read mapping: $name mismatch string contains unexpected symbols\n";
			}
			$string =~s/\^[A-Za-z]+//g; $string =~s/[0-9]+//g; my $num_mis = length($string); # the number of mismatch nucleatides (SNPs+Indels) in read mapping

			next if ( exists($multi_read{$name}) ); # remove the multiple mapped reads
			if ( $one_match eq "scaffold" ) {
				if ( $num_mis < 4 ) { # the default setting -- maximum allowed mismatching number
					push @{$discordant_split{$name}}, [$one_pos, "discordant_split", $one_match, $seq, $other_match, $flag, $quality]; # read across breakpoint
				} else {
					push @{$discordant_split{$name}}, [$one_pos, "mismatch_fail", $one_match, $seq, $other_match, $flag, $quality]; # read with many mismatch base
				}
			} elsif ( $other_match eq "scaffold" ) {
				if ( $num_mis < 4 ) { # the default setting -- maximum allowed mismatching number
					push @{$discordant_split{$name}}, [$other_pos, "discordant_split", $one_match, $seq, $other_match, $flag, $quality]; # read across breakpoint
				} else {
					push @{$discordant_split{$name}}, [$other_pos, "mismatch_fail", $one_match, $seq, $other_match, $flag, $quality]; # read with many mismatch base
				}
			} else {
				print "Step2-2 discoardant-split reads mapping: no hit matches to scaffold (one mapping to GeneA/GeneB; the other mapping to scaffold)\n";
			}
		}
		close IN;

		foreach my $name ( keys %discordant_split ) {
			if ( scalar(@{$discordant_split{$name}}) == 2 ) { # if one end of read show jumpping between scaffold and GeneA / GeneB (in fact -- spanning read);
				if ( $discordant_split{$name}[0][1] eq "discordant_split" and $discordant_split{$name}[1][1] eq "discordant_split" ) {

					# if the $discordant_split{$name}[0] is scaffold sequence
					if ( $discordant_split{$name}[0][2] eq "scaffold" ) {
						my $length_read_t = length($discordant_split{$name}[0][3]); # the true read length -- consider trimming process
						if ( $discordant_split{$name}[1][2] eq $geneA ) { # one-end mapped to upstream gene; the other end mapped to scaffold
							if ( ($discordant_split{$name}[0][0] + $length_read_t) - $breakpoint_scaffold >= $archor ) { # meet anchor requirement 
								if ( ($breakpoint_scaffold - $discordant_split{$name}[0][0]) >= $archor ) { # meet discordant-split read requirement
									$dis_ref->{$name} = $discordant_split{$name}; # correct condition for anchor filtering
								} else {
									$discordant_split{$name}[0][2] = $geneB;
									$dis_ref_span->{$name} = $discordant_split{$name}; # print "rename $name discordant-split reads to spanning reads (geneA ======> <=|===== scaffold(geneB))\n";
								}	
							} else {
								print "Step2-2 discordant split read mapping: $name fails to meet criteria of anchor length (geneA ======> <=====|= scaffold)\n";
							}
						} elsif ( $discordant_split{$name}[1][2] eq $geneB ) { # one-end mapped to downstream gene; the other end mapped to scaffold
							if ( ($breakpoint_scaffold - $discordant_split{$name}[0][0]) >= $archor ) { # meet discordant-split read requirement
								if ( ($discordant_split{$name}[0][0] + $length_read_t) - $breakpoint_scaffold >= $archor ) { # meet anchor requirement 
									$dis_ref->{$name} = $discordant_split{$name}; # correct condition for anchor filtering 
								} else {
									$discordant_split{$name}[0][2] = $geneA;
									$dis_ref_span->{$name} = $discordant_split{$name}; # print "rename $name discordant-split reads to spanning reads (scaffold(geneA)=====|=> <====== geneB)\n";
								}
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria of anchor length: (scaffold =|=====> <====== geneB)\n";
							}
						} else {
							print "Step2-2 discordant-split read mapping: $name matched to geneA or geneB with wrong gene name\n";
						}
					} else { # if $discordant_split{$name}[1] is scaffold sequence
						my $length_read_t = length($discordant_split{$name}[1][3]); # the true read length -- consider trimming process
						if ( $discordant_split{$name}[0][2] eq $geneA ) { # one-end mapped to upstream gene; the other end mapped to scaffold
							if ( ($discordant_split{$name}[1][0] + $length_read_t) - $breakpoint_scaffold >= $archor ) { # meet anchor requirement 
								if ( ($breakpoint_scaffold - $discordant_split{$name}[1][0]) >= $archor ) { # meet discordant-split read requirement
									$dis_ref->{$name} = $discordant_split{$name}; # correct condition for anchor filtering
								} else {
									$discordant_split{$name}[1][2] = $geneB;
									$dis_ref_span->{$name} = $discordant_split{$name}; # print "rename $name discordant-split reads to spanning reads (geneA ======> <=|===== scaffold(geneB))\n";
								}
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria of anchor length: (geneA ======> <=====|= scaffold)\n";
							}
						} elsif ( $discordant_split{$name}[0][2] eq $geneB ) { # one-end mapped to downstream gene; the other end mapped to scaffold
							if ( ($breakpoint_scaffold - $discordant_split{$name}[1][0]) >= $archor ) { # meet discordant-split read requirement
								if ( ($discordant_split{$name}[1][0] + $length_read_t) - $breakpoint_scaffold >= $archor ) { # meet anchor requirement 
									$dis_ref->{$name} = $discordant_split{$name}; # correct condition for anchor filtering 
								} else {
									$discordant_split{$name}[1][2] = $geneA;
									$dis_ref_span->{$name} = $discordant_split{$name}; # print "rename $name discordant-split reads to spanning reads (scaffold(geneA) =====|=> <====== geneB)\n";
								}
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria of anchor length: (scaffold =|=====> <====== geneB)\n";
							}
						} else {
							print "Step2-2 discordant-split read mapping: $name matched to geneA or geneB with wrong gene name\n";
						}
					}
				} elsif ( $discordant_split{$name}[0][1] eq "discordant_split" and $discordant_split{$name}[1][1] eq "mismatch_fail" ) { # only first pass filtering
					if ( $discordant_split{$name}[0][2] eq "scaffold" ) {
						if ( ($discordant_split{$name}[0][0]+$archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($discordant_split{$name}[0][0]+length($discordant_split{$name}[0][3])-$archor) ) {
							$dis_to_sing_ref->{$name} = $discordant_split{$name}; # if one-end mapping with many mismatches, regard it as unmapped.
						} # regard as singleton
					}
				} elsif ( $discordant_split{$name}[1][1] eq "discordant_split" and $discordant_split{$name}[0][1] eq "mismatch_fail" ) { # only second pass filtering
					if ( $discordant_split{$name}[1][2] eq "scaffold" ) {
						if ( ($discordant_split{$name}[1][0]+$archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($discordant_split{$name}[1][0]+length($discordant_split{$name}[1][3])-$archor) ) {
							$dis_to_sing_ref->{$name} = $discordant_split{$name}; # if one-end mapping with many mismatches, regard it as unmapped.
						} # regard as singleton
					}
				} else {
					print "Step2-2 discordant-split read mapping: $name does not meet criteria because of many mismatches in alignment\n";
				}
			}
		}

		my %discordant_split_scaffold; # for both reads mapped to scaffold
		#*****/2. mapping quality >=40; ($7=="=" -- discordant: both reads mapped to scaffold)
		if ( $aligner eq "hisat2" ) {
			open (IN, "awk  -F '\t' -v OFS='\t' '(\$7==\"=\" && \$8>0 && \$5>=40 && \$9!=0)' $sam_file/tmp/noclip.sam | grep 'caffold' |") || die "Step2-2 discordant-split read mapping for hisat2: cannot run awk 2:$!\n";
		} elsif ( $aligner eq "star" ) {
			open (IN, "(awk  -F '\t' -v OFS='\t' '(\$7==\"=\" && \$8>0 && \$5>=3 && \$9!=0)' $sam_file/tmp/noclipAligned.out.sam && awk -F '\t' -v OFS='\t' '(\$3==\$7 && \$5>=3 && \$8>0)' $sam_file/tmp/noclip.sam) | grep 'caffold' |") || die "Step2-2 discordant-split read mapping for star: cannot run awk 2:$!\n";
		}

		while ( <IN> ) {
			chomp $_; #
			my ($name, $flag, $one_match, $one_pos, $other_match, $other_pos, $seq, $quality) = (split /\t/, $_)[0,1,2,3,6,7,9,10];
			$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # trimmed header
			my $string = ""; # record mismatching information
			if ( $_ =~/MD\:Z\:([\w\^]+)/ ) {
				$string = $1;
			} else {
				print "Step2-2 discordant-split read mapping: $name mismatch string contains unexpected symbol\n";
			}
			$string =~s/\^[A-Za-z]+//g; $string =~s/[0-9]+//g; my $num_mis = length($string); # the number of mismatch nucleatides (SNPs+Indels) in read mapping

			if ( $one_match eq "scaffold" ) { 
				if ( $num_mis < 4 ) { # the default setting -- mismatching number
					push @{$discordant_split_scaffold{$name}}, [$one_pos, "discordant_split", $one_match, $seq, $other_match, $flag, $quality]; # read across breakpoint
				} else {
					push @{$discordant_split_scaffold{$name}}, [$one_pos, "mismatch_fail", $one_match, $seq, $other_match, $flag, $quality]; # read with many mismatch base
				}
			} else {
				print "Step2-2 discordant-split read mapping: no hit matches to scaffold (both read mapped to scaffold)\n";
			}
		}
		close IN;


		foreach my $name ( keys %discordant_split_scaffold ) {
			if ( scalar(@{$discordant_split_scaffold{$name}}) == 2 ) {
				if ( $discordant_split_scaffold{$name}[0][1] eq "discordant_split" and $discordant_split_scaffold{$name}[1][1] eq "discordant_split" ) {

					my $length_read_first = length($discordant_split_scaffold{$name}[0][3]); # read length of $discordant_split_scaffold{$name}[0]
					my $length_read_second = length($discordant_split_scaffold{$name}[1][3]); # read length of $discordant_split_scaffold{$name}[1]
					if ( $discordant_split_scaffold{$name}[0][0] < $discordant_split_scaffold{$name}[1][0] ) { # $discordant_split_scaffold{$name}[0] =><= $discordant_split_scaffold{$name}[1]
						my $end_1 = $discordant_split_scaffold{$name}[0][0] + $length_read_first;
						my $end_2 = $discordant_split_scaffold{$name}[1][0] + $length_read_second;

						if ( ($end_1 - $discordant_split_scaffold{$name}[1][0]) >= 4 ) { # if the overlapping region >= 4 bp -- combine two reads as a long fragement
							if ( ($discordant_split_scaffold{$name}[0][0] + $archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($end_2 - $archor) ) {
								$dis_ref->{$name} = $discordant_split_scaffold{$name};
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria of anchor length [0=>1]: (scaffold =|==========|= scaffold)\n";
							}
						} else { # still consider two paired end reads as independent
							my $span_s = $end_1 - $archor; my $span_e = $discordant_split_scaffold{$name}[1][0] + $archor;
							my $split_s = $discordant_split_scaffold{$name}[0][0] + $archor; my $split_e = $end_2 - $archor;
							if ( $split_s <= $breakpoint_scaffold and $breakpoint_scaffold <= $span_s ) { # define as discordant-split read
								$dis_ref->{$name} = $discordant_split_scaffold{$name};
							} elsif ( $span_s < $breakpoint_scaffold and $breakpoint_scaffold < $span_e ) { # define as spanning read
								$discordant_split_scaffold{$name}[0][2] = $geneA;
								$discordant_split_scaffold{$name}[1][2] = $geneB;
								$dis_ref_span->{$name} = $discordant_split_scaffold{$name}; # rename discordant-split reads to spanning reads (scaffold ======||====== geneB)
							} elsif ( $span_e <= $breakpoint_scaffold and $breakpoint_scaffold <= $split_e ) { # define as discordant-split read
								$dis_ref->{$name} = $discordant_split_scaffold{$name};
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria of anchor length [0=>1]: (scaffold =|==========|= scaffold)\n";
							}
						}	
					} else { # $discordant_split_scaffold{$name}[1] =><= $discordant_split_scaffold{$name}[0]
						my $end_1 = $discordant_split_scaffold{$name}[1][0] + $length_read_second;
						my $end_2 = $discordant_split_scaffold{$name}[0][0] + $length_read_first;

						if ( ($end_1 - $discordant_split_scaffold{$name}[0][0]) >= 4 ) { # if the overlapping region >= 4 bp -- combine two reads as a long fragement
							if ( ($discordant_split_scaffold{$name}[1][0] + $archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($end_2 - $archor) ) {
								$dis_ref->{$name} = $discordant_split_scaffold{$name};
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria of anchor length [1=>0]: (scaffold =|==========|= scaffold)\n";
							}
						} else { # still consider two paired end reads as independent
							my $span_s = $end_1 - $archor; my $span_e = $discordant_split_scaffold{$name}[0][0] + $archor;
							my $split_s = $discordant_split_scaffold{$name}[1][0] + $archor; my $split_e = $end_2 - $archor;
							if ( $split_s <= $breakpoint_scaffold and $breakpoint_scaffold <= $span_s ) { # define as discordant-split read
								$dis_ref->{$name} = $discordant_split_scaffold{$name};
							} elsif ( $span_s < $breakpoint_scaffold and $breakpoint_scaffold < $span_e ) { # define as spanning read
								$discordant_split_scaffold{$name}[0][2] = $geneB;
								$discordant_split_scaffold{$name}[1][2] = $geneA;
								$dis_ref_span->{$name} = $discordant_split_scaffold{$name}; # rename discordant-split reads to spanning reads (scaffold ======||====== geneB)
							} elsif ( $span_e <= $breakpoint_scaffold and $breakpoint_scaffold <= $split_e ) { # define as discordant-split read
								$dis_ref->{$name} = $discordant_split_scaffold{$name};
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria of anchor length [1=>0]: (scaffold =|==========|= scaffold)\n";
							}
						}

					}
				} else {
					print "Step2-2 discordant-split read mapping: $name does not meet criteria because of many mismatches in alignment (both scaffold)\n";
				}
			} else {
				if ( $discordant_split_scaffold{$name}[0][1] eq "discordant_split" ) {
					if ( $discordant_split_scaffold{$name}[0][2] eq "scaffold" and $discordant_split_scaffold{$name}[0][4] eq "scaffold" ) {
						my $length_read_first = length($discordant_split_scaffold{$name}[0][3]); #
						if ( ($discordant_split_scaffold{$name}[0][0] + $length_read_first) - $breakpoint_scaffold >= $archor ) { # meet anchor requirement 
							if ( ($breakpoint_scaffold - $discordant_split_scaffold{$name}[0][0]) >= $archor ) { # meet anchor requirement
								if ( $discordant_split_scaffold{$name}[0][5] == 99 || $discordant_split_scaffold{$name}[0][5] == 147 || $discordant_split_scaffold{$name}[0][5] == 83 || $discordant_split_scaffold{$name}[0][5] == 163 ) { # make sure the correct mapping flag
									$dis_ref->{$name} = $discordant_split_scaffold{$name};
								} else {
									print "Step2-2 discordant-split read mapping: $name fails to meet criteria of flag (wrong, only scaffold)\n";
								}
							} else {
								print "Step2-2 discordant-split read mapping: $name fails to meet criteria because of unusual mapping (only scaffold)\n";
							}
						} else {
							print "Step2-2 discordant-split read mapping: $name fails to meet criteria because of unusual mapping (only scaffold)\n";
						}
					} else {
						print "Step2-2 discordant-split read mapping: $name fails to meet criteria because of no scaffold mapping (only scaffold)\n";
					}
				} else {
					print "Step2-2 discordant-split read mapping: $name fails to meet criteria because of many mismatches in alignment (only scaffold)\n";
				}
			}
		}
	}

	sub Singlton {
		my ($dis_ref, $archor, $breakpoint_scaffold, $read_length, $sam_file, $aligner) = @_;
		
		my %spanning_read; # collect the mapped spanning reads for control of singleton read
		# select mapped singlton split read (one mapping to scaffold; the other no mapping hit)
		if ( $aligner eq "hisat2" ) {
			#*****/ mapping quality >= 40; ($9=0 -- singlton mapped; $2=73|89|137|153 -- singlton mapping feature filtering) / for hisat2
			open (IN, "awk -F '\t' -v OFS='\t' '(\$7==\"=\" && \$3==\"scaffold\" && \$5>=40 && \$8>0 && \$9==0 && (\$2==73 || \$2==89 || \$2==137 || \$2==153))' $sam_file/tmp/noclip.sam |") || die "cannot run awk singlton script for hisat2:$!\n";
		} elsif ( $aligner eq "star" ) {
			#*****/ mapping quality >= 0; ($3!=$7 && $7!="=" -- all discordant and spanning reads) / remove singlton read that present in %spanning_read for star
			open (IN, "(awk -F '\t' -v OFS='\t' '(\$3!=\$7 && \$7!=\"=\" && \$5>=0 && \$8>0)' $sam_file/tmp/noclip.sam && awk -F '\t' -v OFS='\t' '(\$3==\$7 && \$5>=3 && \$8>0 && \$3==\"scaffold\")' $sam_file/tmp/noclip.sam) |") || die "Step2-2 spanning read mapping for star: cannot run awk 1:$!\n";
			while ( <IN> ) {
				chomp $_; #
				my ($name, $flag, $one_match, $one_pos, $value, $other_match, $other_pos, $seq, $quality) = (split /\t/, $_)[0,1,2,3,4,6,7,9,10];
				$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # trimmed header
				push @{$spanning_read{$name}}, [$one_pos, "spanning", $one_match, $seq, $other_match, $flag, $quality];
			}
			close IN;
			#*****/ mapping quality > 4; ($9=0 -- singlton mapped; $2=73|89|137|153 -- possible singlton mapping feature filtering) / for start
			open (IN, "awk -F '\t' -v OFS='\t' '(\$7!=\"=\" && \$3==\"scaffold\" && \$5>4 && \$9==0 && (\$2==73 || \$2==89 || \$2==137 || \$2==153))' $sam_file/tmp/noclipAligned.out.sam |") || die "cannot run awk singlton script for star:$!\n";
		}

		while ( <IN> ) {
			chomp $_; #
			my ($name, $flag, $one_match, $one_pos, $other_match, $other_pos, $seq, $quality) = (split /\t/, $_)[0,1,2,3,6,7,9,10];
			$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # trimmed header
			my $string = ""; # record mismatching information
			if ( $_ =~/MD\:Z\:([\w\^]+)/ ) {
				$string = $1;
			} else {
				print "Step2-2 singleton-split read mapping: $name mismatch string contains unexpected symbols\n";
			}
			$string =~s/\^[A-Za-z]+//g; $string =~s/[0-9]+//g; my $num_mis = length($string); # the number of mismatch nucleatides (SNPs+Indels) in read mapping

			next if ( exists($spanning_read{$name}) ); # remove the read that also present in discordant / spanning records
			if ( $one_match eq "scaffold" ) {
				if ( ($one_pos+$archor) <= $breakpoint_scaffold and $breakpoint_scaffold <= ($one_pos+length($seq)-$archor) ) { # read cross breakpoint: get extension (2, -2) -- not applied yet
					if ( length($seq) > 90 and $num_mis < 4 ) { # if mapping region > 90 bp, allow maximum 3 mismatches for singleton
						push @{$dis_ref->{$name}}, [$one_pos, "singlton_split", $one_match, $seq, $flag, $quality];
						# print "singlton: $name\t$singlton{$name}[0][0]\t$singlton{$name}[0][1]\n";
					} elsif ( length($seq) <= 90 and $num_mis < 3 ) { # if mapping region < 90 bp, allow maximum 2 mismatches for singleton
						push @{$dis_ref->{$name}}, [$one_pos, "singlton_split", $one_match, $seq, $flag, $quality];
						# print "singlton: $name\t$singlton{$name}[0][0]\t$singlton{$name}[0][1]\n";
					} else {
						print "Step2-2 singleton-split read mapping: $name fails to meet criteria because of many mismatches in alignment\n";
					}
				} else {
					print "Step2-2 singleton-split read mapping: $name fails to meet criteria because of anchor length: (scaffold =|====|=)\n";
				}
			} else {
				print "Step2-2 singleton-split read mapping: no hit matches to scaffold (one to scaffold; the other unmapped)\n";
			}
		}
		close IN;
	}
1;
