package Check;
use strict;
use warnings;

	sub read_length {
		my ($fastq_1, $fastq_2) = @_;
		my $tmp_1; my $tmp_2; my $seq_1; my $seq_2; # set judgement
		if ( $fastq_1 =~/\.gz/ ) { # compresssed format for the first end
			$tmp_1 = `zcat $fastq_1 | head -n 2`; chomp $tmp_1; $seq_1 = (split /\n/, $tmp_1)[1];
		} else { # uncompressed format for the first end
			$tmp_1 = `head -n 2 $fastq_1`; chomp $tmp_1; $seq_1 = (split /\n/, $tmp_1)[1];
		}
		if ( $fastq_2 =~/\.gz/ ) { # compresssed format for the second end
			$tmp_2 = `zcat $fastq_2 | head -n 2`; chomp $tmp_2; $seq_2 = (split /\n/, $tmp_2)[1];
		} else { # uncompressed format for the second end
			$tmp_2 = `head -n 2 $fastq_2`; chomp $tmp_2; $seq_2 = (split /\n/, $tmp_2)[1];
		}
		my $len_1 = length($seq_1); my $len_2 = length($seq_2);

		if ( $len_1 == $len_2 ) {
			return($len_1);
		} else {
			print "Step 0-1: Read length of fastq file is error\n"; exit;
			return(0);
		}
	}

	sub judge_header {
		my ($fastq) = @_;
		my $tmp; my $name; # set judgement
		if ( $fastq =~/\.gz/ ) { # compresssed format
			$tmp = `zcat $fastq | head -n 2`; chomp $tmp; $name = (split /\n/, $tmp)[0]; # header of first_end
		} else { # uncompressed format
			$tmp = `head -n 2 $fastq`; chomp $tmp; $name = (split /\n/, $tmp)[0]; # header of first_end
		}

		if ( $name =~/\s/ ) {
			return(" ");
		} else {
			if ( $name =~/\// ) {
				return("/");
			} else {
				return("end");
			}
		}
	}

	sub extract_scaffold { # Load scaffold seq; at current stage, only one breakpoint sequence was allowed to involve in alignment 
		my ($read_l, $scaffold, $ref) = @_;
		$/ = ">";
		open (IN, "$scaffold") || die "Step 0-3: Cannot load scaffold sequence dataset:$!\n";
		while ( <IN> ) {
			chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
			my @array = (split /\n/, $_);
			my $header = shift @array; $header =~s/\>//g;
			my $seq = join '', @array;
			my $front; my $back;
			if ( $seq =~/\|/ ) { 
				($front, $back) = (split /\|/, $seq)[0, 1];
			} elsif ( $seq =~/\*/ ) {
				($front, $back) = (split /\*/, $seq)[0, 1];
			} else {
				print "Step 0-3: Breakpoint sequence has wrong split symbol, please use '|' or '*' as breakpoint symbol\n"; exit;
			}
			
			# Compare the length of given breakpoint seq with read length
			my $front_frag; my $front_tag = 0;
			if ( length($front) >= $read_l ) { # if given length of front seq >= read_length; trim length of front seq identical to read_length
				$front_frag = substr($front, -$read_l);
			} else { # if given length of front seq < read_length; mark $front_tag as value of "$read_l - length($front)". 
				$front_frag = $front;
				$front_tag = $read_l - length($front);
			}
			my $back_frag; my $back_tag = 0;
			if ( length($back) >= $read_l ) { #if given length of back seq >= read_length; trim length of back seq identical to read_length
				$back_frag = substr($back, 0, $read_l);
			} else { #if given length of back seq < read_length; mark $back_tag as value of "$read_l - length($back)". 
				$back_frag = $back;
				$back_tag = $read_l - length($back);
			}
			if ( exists($ref->{$header}) ) { print "\nStep 0-3: $header have replicates in scaffold sequence list\n"; }
			$ref->{$header}[0] = [$front_tag, $front_frag]; # breakpoint seq of GeneA
			$ref->{$header}[1] = [$back_tag, $back_frag]; # breakpoint seq of GeneB

			#set reverse complementary sequence
			$front_frag =~tr/ATCG/TAGC/; $front_frag =~tr/[a-z]/[A-Z]/; $front_frag=reverse($front_frag);
			$back_frag =~tr/ATCG/TAGC/; $back_frag =~tr/[a-z]/[A-Z]/; $back_frag=reverse($back_frag);
			$ref->{$header}[2] = [$front_tag, $front_frag]; # reverse complementary breakpoint seq of GeneA
			$ref->{$header}[3] = [$back_tag, $back_frag]; # reverse complementary breakpoint seq of GeneB
			#print "original: $ref->{$header}[0][1] $ref->{$header}[1][1]; reverse_complementary: $ref->{$header}[2][1] $ref->{$header}[3][1]\n"; #--testing
		}
		close IN; 
	}

	sub judge_scaffold {
		my ($scaffold, $tran_seq, $gene, $geneA, $geneB, $resource, $user_path) = @_;

		# Load gene symbol file -- download from ensembl Biomart, Dec 4th 2016
		my %name; 
		my %chrom_include = ("1"=>1,"2"=>1,"3"=>1,"4"=>1,"5"=>1,"6"=>1,"7"=>1,"8"=>1,"9"=>1,"10"=>1,"11"=>1,"12"=>1,"13"=>1,"14"=>1,"15"=>1,"16"=>1,"17"=>1,"18"=>1,"19"=>1,"20"=>1,"21"=>1,"22"=>1,"X"=>1,"Y"=>1,"MT"=>1); # only consider the genes in chromosme (1..22, X, Y and MT); gene within scaffold is filtered out.
		$/ = "\n";
		open (IN, "cut -s -f1,5,10 $gene | sort | uniq |") || die "Step 1: Cannot open gene annotation file:$!\n";
		while ( <IN> ) {
			chomp $_; my ($ensembl, $symbol, $chrom) = (split /\t/, $_)[0, 1, 2];
			if ( exists($chrom_include{$chrom}) ) {
				$name{$symbol} = $ensembl; # key => gene_symbol; value => ensembl_id
			}
		}
		close IN;

		# Load transcript (cdna seq) -- download from ensembl Biomart, Feb 2017
		my %cdna;
		$/ = ">"; # separator symbol
		if ( $resource eq "ensembl" ) {
			open (IN, "$tran_seq") || die "Step 1: Cannot load ensembl transcript sequence dataset:$!\n";
			while ( <IN> ) {
				chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
				my @array = (split /\n/, $_);
				my $header = shift @array; $header =~s/\>//g;
				my ($gene, $tran) = (split /\|/, $header);
				if (! $gene =~/ENSG/ ) { print "Step 1: gene id $gene is wrong\n"; exit; }
				if (! $tran =~/ENST/ ) { print "Step 1: transcript id $tran is wrong\n"; exit; }
				my $seq = join '', @array; my $len = length($seq);
				$cdna{$gene}{$len} = [$tran, $seq]; #print "$gene $len $tran\n$seq\n"; #--testing
			}
			close IN;
		} elsif ( $resource eq "gencode" ) {
			open (IN, "$tran_seq") || die "Step 1: Cannot load gencode transcript sequence dataset:$!\n";
			while ( <IN> ) {
				chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
				my @array = (split /\n/, $_);
				my $header = shift @array; $header =~s/\>//g;
				my ($tran, $gene) = (split /\|/, $header)[0, 1];
				$tran =~s/\.[\d]+$//g; $gene =~s/\.[\d]+$//g;
				if (! $gene =~/ENSG/ ) { print "Step 1: gene id $gene is wrong\n"; exit; }
				if (! $tran =~/ENST/ ) { print "Step 1: transcript id $tran is wrong\n"; exit; }
				my $seq = join '', @array; my $len = length($seq);
				$cdna{$gene}{$len} = [$tran, $seq]; #print "$gene $len $tran\n$seq\n"; #--testing
			}
			close IN;
		} else {
			my %name_tmp; # intermiddle file
			open (UCSC, "cut -s -f2,13 $resource |") || die "Step 1: Cannot load ucsc transcript annotation:$!\n";
			while ( <UCSC> ) {
				chomp $_; my ($tran, $gene) = (split /\t/, $_)[0, 1];
				if ( exists($name{$gene}) ) {
					$name_tmp{$tran} = $name{$gene}; # $trans = transcipt_id; $name{$gene} = ensembl_gene_id
				}
			}
			close UCSC;

			open (IN, "$tran_seq") || die "Step 1: Cannot load ucsc transcript sequence dataset:$!\n";
			while ( <IN> ) {
				chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
				my @array = (split /\n/, $_);
				my $header = shift @array; $header =~s/\>//g;
				my $tran = (split /\s/, $header)[0];
				if ( exists($name_tmp{$tran}) ) {
					my $seq = join '', @array; my $len = length($seq);
					$seq =~tr/atcg/ATCG/;
					$cdna{$name_tmp{$tran}}{$len} = [$tran, $seq]; #print "$gene $len $tran\n$seq\n"; #--testing
				}
			}
			close IN;
		}	

		# see whether user-defined transcript reference sequence present
		if ( $user_path ) {
			open (IN, "$user_path") || die "Step 1: Cannot load user-defined transcript reference sequences:$!\n";
			while ( <IN> ) {
				chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
				my @array = (split /\n/, $_);
				my $header = shift @array; $header =~s/\>//g;
				my ($gene, $tran) = (split /\|/, $header)[0, 1];
				if (! defined($gene) ) { print "Header of user-defined transcript references are error, please follow >gene|transcript\n"; exit; }
				if (! defined($tran) ) { print "Header of user-defined transcript references are error, please follow >gene|transcript\n"; exit; }
				my $seq = join '', @array; my $len = length($seq);
				$cdna{$gene}{$len} = [$tran, $seq];
			}
			close IN;
		} else {
			print "Step1: User-defind transcript reference sequence not present, please continue\n";
		}

		$/ = "\n";
		# Judge whether input gene names are validated in the database; and whether breakpoint sequences match to the cDNA sequences
		my %scaff_final; # %scaff_final data structure: 
				 # $scaff_final{"GeneA"} = [sequence, transcript_id, break_seq_geneB]
				 # $scaff_final{"GeneB"} = [sequence, transcript_id, break_seq_geneB]
		my $tag_scaff_A = 0; # if == 0, scaffold sequence of geneA fails to match cDNA of geneA; if == 1, pass; 
		my $tag_scaff_B = 0; # if == 0, scaffold sequence of geneB fails to match cDNA of geneB; if == 1, pass; 
		if ( $geneA =~/ENSG/ ) { # geneA names as ensembl id
			if ( $geneB =~/ENSG/ )  { # geneB names as ensembl id
				if ( exists($cdna{$geneA}) ) {
					print "Step 1: $geneA is valid Ensembl gene\n";	my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
					$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
				} else {
					print "Step 1: $geneA is not in Ensembl database, please input the valid gene or use user-defined sequences\n";
				}

				if ( exists($cdna{$geneB}) ) {
					print "Step 1: $geneB is valid Ensembl gene\n";	my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
					$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
				} else {
					print "Step 1: $geneB is not in Ensembl database, please input the valid gene or use user-defined sequences\n";
				}
			} else { # geneB names as gene symbol
				if (exists($cdna{$geneA}) ) {
					print "Step 1: $geneA is valid Ensembl gene\n"; my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
					$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
				} else {
					print "Step 1: $geneA is not in Ensembl database, please input the valid gene or use user-defined sequences\n";
				}

				if ( exists($name{$geneB}) ) {
					if ( exists($cdna{$name{$geneB}}) ) {
						print "Step 1: $geneB is valid gene symbol\n"; my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
						$tag_scaff_B = &output_geneB($cdna{$name{$geneB}}, $scaffold, $scaff_final{$geneB});
						if ( $tag_scaff_B != 1 ) {
							if ( exists($cdna{$geneB}) ) {
								$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB}); # define the sub- data structure
							}
						}
					} else {
						if ( exists($cdna{$geneB}) ) {
							print "Step 1: $geneB is user-defined name\n"; my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
							$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
						} else {
							print "Step 1: $geneB is not in gene symbol database, please input the valid gene or use user-defined sequences\n";
						}
					}
				} else {
					if ( exists($cdna{$geneB}) ) {
						print "Step 1: $geneB is user-defined name\n"; my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
						$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
					} else {
						print "Step 1: $geneB is not in gene symbol database, please input the valid gene or use user-defined sequences\n";
					}
				}
			}
		} else { # geneA names as gene symbol
			if ( $geneB =~/ENSG/ )  { # geneB names as ensembl id
				if ( exists($name{$geneA}) ) {
					if ( exists($cdna{$name{$geneA}}) ) {
						print "Step 1: $geneA is valid gene symbol\n"; my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
						$tag_scaff_A = &output_geneA($cdna{$name{$geneA}}, $scaffold, $scaff_final{$geneA});
						if ( $tag_scaff_A != 1 ) {
							if ( exists($cdna{$geneA}) ) {
								$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
							}
						}
					} else {
						if ( exists($cdna{$geneA}) ) {
							print "Step 1: $geneA is user-defined name\n"; my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
							$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
						} else {
							print "Step 1: $geneA is not in gene symbol database, please input the valid gene or use user-defined sequences\n";
						}
					}
				} else {
					if ( exists($cdna{$geneA}) ) {
						print "Step 1: $geneA is user-defined name\n"; my @array; $scaff_final{$geneA} = \@array;
						$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
					} else {
						print "Step 1: $geneA is not in gene symbol database, please input the valid gene or use user-defined sequences\n";
					}
				}

				if ( exists($cdna{$geneB}) ) {
					print "Step 1: $geneB is valid Ensembl gene\n"; my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
					$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
				} else {
					print "Step 1: $geneB is not in Ensembl database, please input the valid gene or use user-defined sequences\n";
				}
			} else { # geneB names as gene symbol
				if ( exists($name{$geneA}) ) {
					if ( exists($cdna{$name{$geneA}}) ) {
						print "Step 1: $geneA is valid gene symbol\n"; my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
						$tag_scaff_A = &output_geneA($cdna{$name{$geneA}}, $scaffold, $scaff_final{$geneA});
						if ( $tag_scaff_A != 1 ) {
							if ( exists($cdna{$geneA}) ) {
								$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
							}
						}
					} else {
						if ( exists($cdna{$geneA}) ) {
							print "Step 1: $geneA is user-defined name\n"; my @array; $scaff_final{$geneA} = \@array;
							$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
						} else {
							print "Step 1: $geneA is not in gene symbol database, please input the valid gene or use user-defined sequences\n";
						}
					}
				} else {
					if ( exists($cdna{$geneA}) ) {
						print "Step 1: $geneA is user-defined name\n"; my @array; $scaff_final{$geneA} = \@array;
						$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
					} else {
						print "Step 1: $geneA is not in gene symbol database, please input the valid gene or use user-defined sequences\n"; 
					}
				}

				if ( exists($name{$geneB}) ) {
					if ( exists($cdna{$name{$geneB}}) ) {
						print "Step 1: $geneB is valid gene symbol\n"; my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
						$tag_scaff_B = &output_geneB($cdna{$name{$geneB}}, $scaffold, $scaff_final{$geneB});
						if ( $tag_scaff_B != 1 ) {
							if ( exists($cdna{$geneB}) ) {
								$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
							}
						}
					} else {
						if ( exists($cdna{$geneB}) ) {
							print "Step 1: $geneB is user-defined name\n"; my @array; $scaff_final{$geneB} = \@array;
							$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
						} else {
							print "Step 1: $geneB is not in gene symbol database, please input the valid gene or use user-defined sequences\n";
						}
					}
				} else {
					if ( exists($cdna{$geneB}) ) {
						print "Step 1: $geneB is user-defined name\n"; my @array; $scaff_final{$geneB} = \@array;
						$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
					} else {
						print "Step 1: $geneB is not in gene symbol database, please input the valid gene or use user-defined sequences\n";
					}
				}
			}
		}

		return($tag_scaff_A, $tag_scaff_B, \%scaff_final);
	}

	sub output_geneB {
		my ($cdna_seq, $scaff_ref, $scaff_seq) = @_;
		
		my $tag = 0; my %hash; # collect all cdna with matched breakpoint seq for geneB 
		foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

			my $split_tag = $scaff_ref->[1][1];
			my @cdna_geneB = split /$split_tag/, $cdna_seq->{$rank}[1]; # print ">$cdna_seq->{$rank}[0]\n$cdna_seq->{$rank}[1]\n";
			if ( scalar(@cdna_geneB) == 2 ) {
				$tag = 1;
				if ( $scaff_ref->[1][0] == 0 ) { # length of breakpoint == length of read
					my $len_geneB = length($scaff_ref->[1][1]);
					push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
				} else {
					my $add = substr($cdna_geneB[1], 0, $scaff_ref->[1][0]); $add = $scaff_ref->[1][1].$add;
					my $len_geneB = length($add);
					push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneB]
				}
				next;
			} elsif ( scalar(@cdna_geneB) >= 3 ) {
				print "Step 1: Breakpoint homolog in fusion sequence of downstream gene partner ($rank)\n";
			} else {
				my $len = length($split_tag) - 3;
				$split_tag = substr($split_tag, -$len); # trim the frist 3 base of breakpoint sequence
				@cdna_geneB = split /$split_tag/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneB) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[1][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneB = length($scaff_ref->[1][1]);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
					} else {
						my $add = substr($cdna_geneB[1], 0, $scaff_ref->[1][0]); $add = $scaff_ref->[1][1].$add;
						my $len_geneB = length($add);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneB]
					}
					next;
				} elsif ( scalar(@cdna_geneB) >= 3 ) {
					print "Step 1: Breakpoint homolog in fusion sequence of downstream gene partner ($rank)\n";
				}
			}
		}

		if (! %hash ) {
			$tag = 0; # re-set $tag = 0;
			foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

				my $split_tag_r = $scaff_ref->[3][1];
				my @cdna_geneB_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneB_r) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[3][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneB = length($scaff_ref->[1][1]);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
					} else {
						my $add = substr($cdna_geneB_r[0], -$scaff_ref->[3][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $scaff_ref->[1][1].$add;
						my $len_geneB = length($add);
						push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; #[sequence, transcript_id, break_seq_geneB]
					}
					next;
				} else {
					my $len_r = length($split_tag_r) - 3;
					$split_tag_r = substr($split_tag_r, 0, $len_r); # trim the last 3 base of breakpoint sequence
					@cdna_geneB_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
					if ( scalar(@cdna_geneB_r) == 2 ) {
						$tag = 1;
						if ( $scaff_ref->[3][0] == 0 ) { # length of breakpoint == length of read
							my $len_geneB = length($scaff_ref->[1][1]);
							push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[1][1]]; # [sequence, transcript_id, break_seq_geneB]
						} else {
							my $add = substr($cdna_geneB_r[0], -$scaff_ref->[3][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $scaff_ref->[1][1].$add;
							my $len_geneB = length($add);
							push @{$hash{$len_geneB}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneB]
						}
						next;
					}
				}
			}
		}

		if ( $tag == 1 ) {
			foreach my $len ( reverse sort {$a <=> $b} keys %hash ) {
				foreach my $ref ( @{$hash{$len}} ) { # select the first element [0] of %hash (the longest transcript match to breakpoint region of geneB)
					$scaff_seq->[0] = $ref->[0]; $scaff_seq->[1] = $ref->[1]; $scaff_seq->[2] = $ref->[2]; # [sequence, transcript_id, break_seq_geneB]
					return(1);
				}
			}
		} else {
			return(0);
		}
	}

	sub output_geneA {
		my ($cdna_seq, $scaff_ref, $scaff_seq) = @_;

		my $tag = 0; my %hash; # collect all cdna with matched breakpoint seq for geneA
		foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

			my $split_tag = $scaff_ref->[0][1];
			my @cdna_geneA = split /$split_tag/, $cdna_seq->{$rank}[1]; # print ">$cdna_seq->{$rank}[0]\n$cdna_seq->{$rank}[1]\n";
			if ( scalar(@cdna_geneA) == 2 ) {
				$tag = 1;
				if ( $scaff_ref->[0][0] == 0 ) { # length of breakpoint == length of read
					my $len_geneA = length($scaff_ref->[0][1]);
					push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
				} else { # length of breakpoint < length of read; need to get extension
					my $add = substr($cdna_geneA[0], -$scaff_ref->[0][0]); $add = $add.$scaff_ref->[0][1]; 
					my $len_geneA = length($add);
					push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
				}
				next;
			} elsif ( scalar(@cdna_geneA) >= 3 ) {
				print "Step 1: Breakpoint homolog in fusion sequence of upstream gene partner ($rank)\n";
			} else {
				my $len = length($split_tag) - 3; 
				$split_tag = substr($split_tag, 0, $len); # trim the last 3 base of breakpoint sequence
				@cdna_geneA = split /$split_tag/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneA) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[0][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneA = length($scaff_ref->[0][1]);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
					} else {
						my $add = substr($cdna_geneA[0], -$scaff_ref->[0][0]); $add = $add.$scaff_ref->[0][1];
						my $len_geneA = length($add);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
					}
					next;
				} elsif ( scalar(@cdna_geneA) >= 3 ) {
					print "Step 1: Breakpoint homolog in fusion sequence of upstream gene partner ($rank)\n";
				}
			}
		}

		if (! %hash) {
			$tag = 0;
			foreach my $rank ( reverse sort {$a <=> $b} keys %{$cdna_seq} ) { # loop from the longest transcript to the shortest transcript

				my $split_tag_r = $scaff_ref->[2][1];
				my @cdna_geneA_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
				if ( scalar(@cdna_geneA_r) == 2 ) {
					$tag = 1;
					if ( $scaff_ref->[2][0] == 0 ) { # length of breakpoint == length of read
						my $len_geneA = length($scaff_ref->[0][1]);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
					} else { # length of breakpoint < length of read; need to get extension
						my $add = substr($cdna_geneA_r[1], 0, $scaff_ref->[2][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $add.$scaff_ref->[0][1];
						my $len_geneA = length($add);
						push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
					}
					next;
				} else {
					my $len_r = length($split_tag_r) - 3; 
					$split_tag_r = substr($split_tag_r, -$len_r); # trim the first 3 base of breakpoint sequence
					@cdna_geneA_r = split /$split_tag_r/, $cdna_seq->{$rank}[1];
					if ( scalar(@cdna_geneA_r) == 2 ) {
						$tag = 1;
						if ( $scaff_ref->[2][0] == 0 ) { # length of breakpoint == length of read
							my $len_geneA = length($scaff_ref->[0][1]);
							push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $scaff_ref->[0][1]]; # [sequence, transcript_id, break_seq_geneA]
						} else { # length of breakpoint < length of read; need to get extension
							my $add = substr($cdna_geneA_r[1], 0, $scaff_ref->[2][0]); $add =~tr/ATCG/TAGC/; $add = reverse($add); $add = $add.$scaff_ref->[0][1];
							my $len_geneA = length($add);
							push @{$hash{$len_geneA}}, [$cdna_seq->{$rank}[1], $cdna_seq->{$rank}[0], $add]; # [sequence, transcript_id, break_seq_geneA]
						}
						next;
					}
				}
			}
		}

		if ( $tag == 1 ) {
			foreach my $len ( reverse sort {$a <=> $b} keys %hash ) {
				foreach my $ref ( @{$hash{$len}} ) { # select the first element [0] of %hash (the longest transcript match to breakpoint region of geneA)
					$scaff_seq->[0] = $ref->[0]; $scaff_seq->[1] = $ref->[1]; $scaff_seq->[2] = $ref->[2]; # [sequence, transcript_id, break_seq_geneA]
					return(1);
				}
			}
		} else {
			return(0);
		}
	}
1;
