package Genome_align;
use strict;
use warnings;

	sub discordant_specif {
		my ($dis_ref, $multiple_ref, $path, $gene, $read_type, $sep, $a_pos_tag, $b_pos_tag) = @_; #read_type: spanning / discordant
		foreach my $id ( keys %{$dis_ref} ) {
			my $com = "@".$id."/"; #please do not use $sep, we have to set the header in end with '/'
			my $hit_1 = `grep "$com" -A 3 $path/${read_type}_1.txt`; chomp $hit_1;
			my $hit_2 = `grep "$com" -A 3 $path/${read_type}_2.txt`; chomp $hit_2;
			my $seq_1 = (split /\n/, $hit_1)[1]; my $seq_1_tr = $seq_1; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
			my $seq_2 = (split /\n/, $hit_2)[1]; my $seq_2_tr = $seq_2; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
			$dis_ref->{$id}[0][1] = $seq_1;
			$dis_ref->{$id}[0][2] = $seq_1_tr; #reverse complementary
			$dis_ref->{$id}[0][3] = $hit_1;
			$dis_ref->{$id}[1][1] = $seq_2;
			$dis_ref->{$id}[1][2] = $seq_2_tr; #reverse complementary
			$dis_ref->{$id}[1][3] = $hit_2;
		}

		# collect mutiple hit reads or bad quality reads
		open (IN, "awk -F '\t' -v OFS='\t' '(\$8>=0 && \$5<40)' $path/tmp/${read_type}_sec.sam |") || die "Step 3: cannot $read_type sam file:$!\n";
		while ( <IN> ) {
			chomp $_; my ($id, $flag, $chr, $pos, $quality, $cigar, $seq) = (split /\t/, $_)[0,1,2,3,4,5,9];
			$id =~s/\/[\w\:\-]+$//g; $id =~s/\s[\w\:\-]+$//g; # print "#*$id\n"; # trimmed header
			push @{$multiple_ref->{$id}{$seq}}, [$chr, $pos, $quality, $cigar, $flag]; #
		}
		close IN;
		# filter out the true positives in the mutiple hit reads
		foreach my $id ( keys %{$multiple_ref} ) {
			my @tmp_n = keys %{$multiple_ref->{$id}};
			if ( scalar(@tmp_n) == 1 ) {
				my $tmp_pos_start; my $tmp_pos_end; my $j = 0; my $tmp_tag = 1;
				for (my $i = 0; $i < scalar(@{$multiple_ref->{$id}{$tmp_n[0]}}); $i++ ) {
					my $tmp_cigar = $multiple_ref->{$id}{$tmp_n[0]}[$i][3]; # get cigar value

					if ( $j == 0 ) { 
						$tmp_pos_start = $multiple_ref->{$id}{$tmp_n[0]}[$i][1];	$tmp_pos_end = $tmp_pos_start;
						if ( $tmp_cigar =~/^([\d]+)S/ ) { # front soft-clip
							if ( $1 > 9 ) { next; } # if soft-clip too many (>9), not trusted for mapping specicifity.
						}
						if ( $tmp_cigar =~/([\d]+)S$/ ) { # back soft-clip
							if ( $1 > 9 ) { # if soft-clip too many (>9), not trusted for mapping specicifity.
								next; 
							} else { 
								$tmp_pos_end = $tmp_pos_end + $1; 
							}
						}
						while ( $tmp_cigar =~/([\d]+)M/gs ) { $tmp_pos_end = $tmp_pos_end + $1; }
						while ( $tmp_cigar =~/([\d]+)N/gs ) { $tmp_pos_end = $tmp_pos_end + $1; }
						while ( $tmp_cigar =~/([\d]+)D/gs ) { $tmp_pos_end = $tmp_pos_end + $1; }
						$j = 1;
					} else {
						my $tmp_comp_start = $multiple_ref->{$id}{$tmp_n[0]}[$i][1];	my $tmp_comp_end = $tmp_comp_start;
						if ( $tmp_cigar =~/^([\d]+)S/ ) { # front soft-clip
							if ( $1 > 9 ) { next; } # if soft-clip too many (>9), not trusted for mapping specicifity.
						}
						if ( $tmp_cigar =~/([\d]+)S$/ ) { # back soft-clip
							if ( $1 > 9 ) { # if soft-clip too many (>9), not trusted for mapping specicifity.
								next; 
							} else {
								$tmp_comp_end = $tmp_comp_end + $1; 
							}
						} 
						while ( $tmp_cigar =~/([\d]+)M/gs ) { $tmp_comp_end = $tmp_comp_end + $1; }
						while ( $tmp_cigar =~/([\d]+)N/gs ) { $tmp_comp_end = $tmp_comp_end + $1; }
						while ( $tmp_cigar =~/([\d]+)D/gs ) { $tmp_comp_end = $tmp_comp_end + $1; }

						if ( $tmp_pos_start != $tmp_comp_start ) {
							if ( $tmp_pos_end != $tmp_comp_end ) {
								$tmp_tag = 0; last; # multiple mapping positive 
							}
						}
					}
				}
				if ( $tmp_tag == 1 ) { # not real multiple mapping 
					delete($multiple_ref->{$id}{$tmp_n[0]});
					delete($multiple_ref->{$id});
				}
			}
		}

		my %mapped; # mapped read with good quality (and remove the other end of reads with multiple hits)
		open (IN, "awk -F '\t' -v OFS='\t' '(\$8>=0 && \$5>=40)' $path/tmp/${read_type}_sec.sam |") || die "Step 3: cannot $read_type sam file:$!\n";
		while ( <IN> ) {
			chomp $_; my ($id, $flag, $chr, $pos, $quality, $cigar, $seq) = (split /\t/, $_)[0,1,2,3,4,5,9];
			$id =~s/\/[\w\:\-]+$//g; $id =~s/\s[\w\:\-]+$//g; # print "#*$id\n"; # trimmed header
			next if ( exists($multiple_ref->{$id}) );
			push @{$mapped{$id}{$seq}}, [$chr, $pos, $quality, $cigar, $flag]; #
		}
		close IN;

		foreach my $id ( keys %{$dis_ref} ) {
			if ( exists($mapped{$id}) ) {
				my $tag = 1; # judgement: 1 => fail; 0 => accept
				my @seq_num = keys %{$mapped{$id}}; 
				if ( scalar(@seq_num) > 2 ) { delete($dis_ref->{$id}); next; } # more than two types of mapping sequence at genomic level
 
				foreach my $seq ( keys %{$mapped{$id}} ) {
					# read with multiple mapping positions (>=2)
					if ( scalar(@{$mapped{$id}{$seq}}) > 1 ) {
						if ( $mapped{$id}{$seq}[0][2] <= $mapped{$id}{$seq}[1][2] ) { # read has multiple hits with equal quality - filtered out

							my $tmp_pos_start; my $tmp_pos_end; my $tmp_qual; my $j = 0; my $tmp_tag = 1;
							for ( my $i = 0; $i < scalar(@{$mapped{$id}{$seq}}); $i++ ) {
								my $tmp_cigar = $mapped{$id}{$seq}[$i][3]; # get cigar value

								if ( $j == 0 ) { 
									$tmp_pos_start = $mapped{$id}{$seq}[$i][1];	$tmp_pos_end = $tmp_pos_start;
									if ( $tmp_cigar =~/^([\d]+)S/ ) { # front soft-clip
										if ( $1 > 9 ) { next; } # if soft-clip too many (>9), not trusted for mapping specicifity.
									}
									if ( $tmp_cigar =~/([\d]+)S$/ ) { # back soft-clip
										if ( $1 > 9 ) {	# soft-clip too many (>9), not trusted for specicifity.
											next; 
										} else { 
											$tmp_pos_end = $tmp_pos_end + $1; 
										}
									}
									while ( $tmp_cigar =~/([\d]+)M/gs ) { $tmp_pos_end = $tmp_pos_end + $1; }
									while ( $tmp_cigar =~/([\d]+)N/gs ) { $tmp_pos_end = $tmp_pos_end + $1; }
									while ( $tmp_cigar =~/([\d]+)D/gs ) { $tmp_pos_end = $tmp_pos_end + $1; }
									$j = 1;
								} else { 
									my $tmp_comp_start = $mapped{$id}{$seq}[$i][1];    my $tmp_comp_end = $tmp_comp_start;
									if ( $tmp_cigar =~/^([\d]+)S/ ) { # front soft-clip
										if ( $1 > 9 ) { next; } # if soft-clip too many (>9), not trusted for mapping specicifity.
									}
									if ( $tmp_cigar =~/([\d]+)S$/ ) { # back soft-clip
										if ( $1 > 9 ) { # soft-clip too many (>9), not trusted for mapping specicifity.
											next;
										} else {
											$tmp_comp_end = $tmp_comp_end + $1; 
										}
									}
									while ( $tmp_cigar =~/([\d]+)M/gs ) { $tmp_comp_end = $tmp_comp_end + $1; }
									while ( $tmp_cigar =~/([\d]+)N/gs ) { $tmp_comp_end = $tmp_comp_end + $1; }
									while ( $tmp_cigar =~/([\d]+)D/gs ) { $tmp_comp_end = $tmp_comp_end + $1; }
	
									if ( $tmp_pos_start != $tmp_comp_start ) {
										if ( $tmp_pos_end != $tmp_comp_end ) {
											$tmp_tag = 0; last; # multiple mapping positive
										}
									}
								}
							}
							if ( $tmp_tag == 0 ) { #multiple mapping positive
								delete($dis_ref->{$id}); last;
							}
						}	
					}
					# read not primary sequence mapped - filtered out
					if ( $mapped{$id}{$seq}[0][4] >= 256 ) { 
						#print "Step 3: ${read_type}-filtered: $id mapped not primary sequence: $seq\n";
						delete($dis_ref->{$id}); last;
					}
					# check whether read mapped to GeneA/GeneB if it has unique mapping position (tag == 1, filtered out)
					if ( $a_pos_tag == 1 ) { # geneA has gene structure annotation
						if ( $b_pos_tag == 1 ) { # geneB has gene structure annotation
							if ( $seq eq $dis_ref->{$id}[0][1] ) {
								&discordant_func_two($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "Two $read_type: $id - $dis_ref->{$id}[0][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[0][2] ) {
								&discordant_func_two($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "Two $read_type: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
							}
							if ( $seq eq $dis_ref->{$id}[1][1] ) {
								&discordant_func_two($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "Two $read_type: $id - $dis_ref->{$id}[1][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
								&discordant_func_two($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "Two $read_type: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
							}
						} else { # Only geneA has gene structure annotation
							if ( $seq eq $dis_ref->{$id}[0][1] ) {
								&discordant_func_one($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - $dis_ref->{$id}[0][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[0][2] ) {
								&discordant_func_one($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
							}
							if ( $seq eq $dis_ref->{$id}[1][1] ) {
								&discordant_func_one($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - $dis_ref->{$id}[1][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
								&discordant_func_one($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
							}
						}
					} else {
						if ( $b_pos_tag == 1 ) { # Only geneB has gene structure annotation
							if ( $seq eq $dis_ref->{$id}[0][1] ) {
								&discordant_func_one($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - $dis_ref->{$id}[0][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[0][2] ) {
								&discordant_func_one($id, $dis_ref->{$id}[0][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
							}
							if ( $seq eq $dis_ref->{$id}[1][1] ) {
								&discordant_func_one($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - $dis_ref->{$id}[1][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
								&discordant_func_one($id, $dis_ref->{$id}[1][0], $mapped{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One $read_type: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
							}
						} else { # this condition is for user-defined transcript sequence input
							$tag = 0;
							#print "Step 3: ${read_type} - $id ($dis_ref->{$id}[0][0]-$dis_ref->{$id}[1][0]) show no anntation in the database, but unique mapping and probably correct\n";
						}
					}

					if ( $tag == 1 ) {
						delete($dis_ref->{$id}); last;
					}
				}
			} else {
				if ( exists($multiple_ref->{$id}) ) {
					#print "Step 3: ${read_type}-filtered: $id shows unspecific mapping in genome alignment\n";
					delete($dis_ref->{$id});
				} else {
					#print "Step 3: ${read_type} - $id ($dis_ref->{$id}[0][0]-$dis_ref->{$id}[1][0]) show no mapping hits in genome alignment, but probably correct\n"; # no-splicing alignment is more sensitive than splicing
				}
			}
		}
	}

	sub singlton_specif {
		my ($dis_ref, $multiple_ref, $path, $gene, $read_type, $sep, $a_pos_tag, $b_pos_tag) = @_;
		foreach my $id ( keys %{$dis_ref} ) {
			my $com = "@".$id."/"; # Please do not use $sep, we have to set the header in end with '/'
			my $hit_1 = `grep "$com" -A 3 $path/${read_type}_1.txt`; chomp $hit_1;
			my $hit_2 = `grep "$com" -A 3 $path/${read_type}_2.txt`; chomp $hit_2;
			my $seq_1 = (split /\n/, $hit_1)[1]; my $seq_1_tr = $seq_1; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
			my $seq_2 = (split /\n/, $hit_2)[1]; my $seq_2_tr = $seq_2; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
			$dis_ref->{$id}[0][1] = $seq_1;
			$dis_ref->{$id}[0][2] = $seq_1_tr; #reverse complementary
			$dis_ref->{$id}[0][3] = $hit_1;
			$dis_ref->{$id}[1][1] = $seq_2;
			$dis_ref->{$id}[1][2] = $seq_2_tr; #reverse complementary
			$dis_ref->{$id}[1][3] = $hit_2;
		}
		# collect all the singlton split reads
		open (IN, "awk  -F '\t' -v OFS='\t' '(\$8>=0)' $path/tmp/${read_type}_sec.sam |") || die "cannot ${read_type} sam file:$!\n";
		while ( <IN> ) {
			chomp $_; my ($name, $flag, $chr, $pos, $quality, $cigar, $seq) = (split /\t/, $_)[0,1,2,3,4,5,9];
			$name =~s/\/[\w\:\-]+$//g; $name =~s/\s[\w\:\-]+$//g; # print "#*$name\n"; # trimmed header
			push @{$multiple_ref->{$name}{$seq}}, [$chr, $pos, $quality, $cigar, $flag];
		}
		close IN;

		foreach my $id ( keys %{$dis_ref} ) {
			if ( exists($multiple_ref->{$id}) ) {
				my $tag = 1; # judgememnt: 1 => fail; 0 => accept
				foreach my $seq ( keys %{$multiple_ref->{$id}} ) {
					# read with multiple mapping positions (>=2)
					if ( scalar(@{$multiple_ref->{$id}{$seq}}) > 1 ) {
						if ( $multiple_ref->{$id}{$seq}[0][2] <= $multiple_ref->{$id}{$seq}[1][2] ) { # read has multiple hits with equal quality - filtered out
							#print "Step 3: ${read_type}-filtered: $id show multiple hits with equal quality: $seq\n"; 
							delete($dis_ref->{$id}); last;
						}
					}
					# read not primary sequence mapped - filtered out
					if ( $multiple_ref->{$id}{$seq}[0][4] >= 256 ) { 
						#print "Step 3: ${read_type}-filtered: $id mapped not primary sequence: $seq\n";
						delete($dis_ref->{$id}); last;
					}
					
					if ( $a_pos_tag == 1 ) {
						if ( $b_pos_tag == 1 ) {
							if ( $seq eq $dis_ref->{$id}[0][1] ) {
								&single_func_two($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "Two singlton: $id - $dis_ref->{$id}[0][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[0][2] ) { 
								my $name_1 = $dis_ref->{$id}[0][0];
								&single_func_two($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
								if ( $name_1 eq "NULL" ) {
									$dis_ref->{$id}[0][0] = $dis_ref->{$id}[0][0]."(reverse_complement)";
								} 
								# print "Two singlton: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
							}
							if ( $seq eq $dis_ref->{$id}[1][1] ) { 
								&single_func_two($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "Two singlton: $id - $dis_ref->{$id}[1][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
								my $name_2 = $dis_ref->{$id}[1][0];
								&single_func_two($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
								if ( $name_2 eq "NULL" ) {
									$dis_ref->{$id}[1][0] = $dis_ref->{$id}[1][0]."(reverse_complement)";
								} 
								# print "Two singlton: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
							}
						} else {
							if ( $seq eq $dis_ref->{$id}[0][1] ) {
								&single_func_one($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One singlton: $id - $dis_ref->{$id}[0][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[0][2] ) {
								my $name_1 = $dis_ref->{$id}[0][0];
								&single_func_one($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
								if ( $name_1 eq "NULL" ) {
									$dis_ref->{$id}[0][0] = $dis_ref->{$id}[0][0]."(reverse_complement)";
								} 
								# print "One singlton: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
							}
							if ( $seq eq $dis_ref->{$id}[1][1] ) { 
								&single_func_one($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One singlton: $id - $dis_ref->{$id}[1][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
								my $name_2 = $dis_ref->{$id}[1][0];
								&single_func_one($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
								if ( $name_2 eq "NULL" ) {
									$dis_ref->{$id}[1][0] = $dis_ref->{$id}[1][0]."(reverse_complement)";
								} 
								# print "One singlton: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
							}
						}
					} else {
						if ( $b_pos_tag == 1 ) {
							if ( $seq eq $dis_ref->{$id}[0][1] ) {
								&single_func_one($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One singlton: $id - $dis_ref->{$id}[0][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[0][2] ) {
								my $name_1 = $dis_ref->{$id}[0][0];
								&single_func_one($id, \$dis_ref->{$id}[0][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
								if ( $name_1 eq "NULL" ) {
									$dis_ref->{$id}[0][0] = $dis_ref->{$id}[0][0]."(reverse_complement)";
								} 
								# print "One singlton: $id - reverse complementary $dis_ref->{$id}[0][0]\t$tag\n";
							}
							if ( $seq eq $dis_ref->{$id}[1][1] ) {
								&single_func_one($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); 
								# print "One singlton: $id - $dis_ref->{$id}[1][0]\t$tag\n";
							} elsif ( $seq eq $dis_ref->{$id}[1][2] ) {
								my $name_2 = $dis_ref->{$id}[1][0];
								&single_func_one($id, \$dis_ref->{$id}[1][0], $multiple_ref->{$id}{$seq}, $gene, \$tag, $read_type); #Note the use the reference for scalar (\$dis_ref->{$id}[0][0]; \$tag)
								if ( $name_2 eq "NULL" ) {
									$dis_ref->{$id}[1][0] = $dis_ref->{$id}[1][0]."(reverse_complement)";
								} 
								# print "One singlton: $id - reverse complementary $dis_ref->{$id}[1][0]\t$tag\n";
							}
						} else {
							$tag = 0;
							#print "Step 3: ${read_type} - $id ($dis_ref->{$id}[0][0]-$dis_ref->{$id}[1][0]) show no anntation in the database, but unique mapping and probably correct\n";
						}
					}
					if ( $tag == 1 ) { delete($dis_ref->{$id}); last; }
				}
			} else {
				#print "Step 3: singlton: $id - not mapped in genome\n";
			}
		}
	}

#####################################
# substitune function within module #
#####################################
	sub discordant_func_one {
		my ($id_ref, $read_ref, $mapped_ref, $partner_ref, $tag_ref, $read_type) = @_;
		my $gene_name = $read_ref; $gene_name = (split /\(/, $gene_name)[0]; $gene_name =~s/\s//g; # edit gene name

		if ( $gene_name eq $partner_ref->[0][3] ) { # if gene_name (A|B) has the gene model annotation
			if ( $mapped_ref->[0][0] eq $partner_ref->[0][0] ) { # the chromosome of mapping read == that of gene
				if ( $partner_ref->[0][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within gene
					$$tag_ref = 0;
				} else {
					#print "Step 3: ${read_type}-filtered: $id_ref ($read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[0][3]\n";
					$$tag_ref = 1;
				}
			}
		} else {
			$$tag_ref = 0;
			#print "Step 3: ${read_type}: $id_ref ($read_ref) has no annotation in database, but unique mapping and probably correct\n";
		}
	}
		
	sub discordant_func_two {
		my ($id_ref, $read_ref, $mapped_ref, $partner_ref, $tag_ref, $read_type) = @_; # Note: $tag_ref is a refrerence

		if ( $partner_ref->[0][0] eq $partner_ref->[1][0] ) { # geneA and geneB are in the same chromosome
			if ( $mapped_ref->[0][0] eq $partner_ref->[0][0] ) {
				if ( $partner_ref->[0][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					$$tag_ref = 0;
				} elsif ( $partner_ref->[1][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					$$tag_ref = 0;
				} else {
					#print "Step 3: ${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
					$$tag_ref = 1;
				}
			} else {
				#print "Step 3: ${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
				$$tag_ref = 1;
			}
		} else { # geneA and geneB are in the different chromosome
			if ( $mapped_ref->[0][0] eq $partner_ref->[0][0] ) { # the chromosome of mapping read == that of geneA
				if ( $partner_ref->[0][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					$$tag_ref = 0;
				} else {
					#print "Step 3: ${read_type}-filtered: $id_ref ($read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[0][3]\n";
					$$tag_ref = 1;
				}
			} elsif ( $mapped_ref->[0][0] eq $partner_ref->[1][0] ) { # the chromosome of mapping read == that of geneB
				if ( $partner_ref->[1][1] < $mapped_ref->[0][1] and $mapped_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					$$tag_ref = 0;
				} else {
					#print "Step 3: ${read_type}-filtered: $id_ref ($read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[1][3]\n";
					$$tag_ref = 1;
				}
			} else {
				#print "Step 3: ${read_type}-filtered: $id_ref ($read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
				$$tag_ref = 1;
			}
		}
	}

	sub single_func_one {
		my ($id_ref, $read_ref, $multiple_ref, $partner_ref, $tag_ref, $read_type) = @_; #Note: $read_ref and $tag_ref are two references

		if ( $multiple_ref->[0][2] < 40 ) { # set mapping quality filtering
			#print "Step 3: $read_type-filtered: $id_ref ($$read_ref) with bad mapping quality\n";
			$$tag_ref = 1; return;
		} 

		if ( $multiple_ref->[0][0] eq $partner_ref->[0][0] ) {
			if ( $partner_ref->[0][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
				if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[0][3]; } # replace geneA to NULL 
				$$tag_ref = 0;
			} else {
				$$tag_ref = 0;
				#print "Step 3: ${read_type}: $id_ref ($read_ref) has no annotation in database, but unique mapping and probably correct\n";
			}
		} else {
			$$tag_ref = 0;
			#print "Step 3: ${read_type}: $id_ref ($read_ref) has no annotation in database, but unique mapping and probably correct\n";
		}
	}	

	sub single_func_two { # whether unmapped read in the frist round (no-splicing) was mapped to geneA or geneB in the second round (splicing)
		my ($id_ref, $read_ref, $multiple_ref, $partner_ref, $tag_ref, $read_type) = @_; #Note: $read_ref and $tag_ref are two references

		if ( $multiple_ref->[0][2] < 40 ) { # set mapping quality filtering
			#print "Step 3: $read_type-filtered: $id_ref ($$read_ref) with bad mapping quality\n";
			$$tag_ref = 1; return;
		}
	
		if ( $partner_ref->[0][0] eq $partner_ref->[1][0] ) { # geneA and geneB are in the same chromosome
			if ( $multiple_ref->[0][0] eq $partner_ref->[0][0] ) {
				if ( $partner_ref->[0][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[0][3]; } # replace geneA to NULL
					$$tag_ref = 0;
				} elsif ( $partner_ref->[1][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[1][3]; } # replace geneB to NULL
					$$tag_ref = 0;
				} else {
					#print "Step 3: ${read_type}-filtered: $id_ref ($$read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
					$$tag_ref = 1;
				}
			} else {
				#print "Step 3: ${read_type}-filtered: $id_ref ($$read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
				$$tag_ref = 1;
			}
		} else {
			if ( $multiple_ref->[0][0] eq $partner_ref->[0][0] ) { # the chromosome of mapping read == that of geneA
				if ( $partner_ref->[0][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[0][2] ) { # mapping position within geneA
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[0][3]; } # replace geneA to NULL
					$$tag_ref = 0;
				} else {
					#print "Step 3: ${read_type}-filtered: $id_ref ($$read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[0][3]\n";
					$$tag_ref = 1;
				}
			} elsif ( $multiple_ref->[0][0] eq $partner_ref->[1][0] ) { #the chromosome of mapping read == that of geneB
				if ( $partner_ref->[1][1] < $multiple_ref->[0][1] and $multiple_ref->[0][1] < $partner_ref->[1][2] ) { # mapping position within geneB
					if ( $$read_ref eq "NULL" ) { $$read_ref = $partner_ref->[1][3]; } # replace geneB to NULL
					$$tag_ref = 0;
				} else {
					#print "Step 3: ${read_type}-filtered: $id_ref ($$read_ref) read mapped the same chrom but wrong gene partners, not $partner_ref->[1][3]\n";
					$$tag_ref = 1;
				}
			} else {
				#print "Step 3: ${read_type}-filtered: $id_ref ($$read_ref) read mapped the wrong gene partners, not $partner_ref->[0][3]/$partner_ref->[1][3]\n";
				$$tag_ref = 1;
			}
		}
	}
1;
