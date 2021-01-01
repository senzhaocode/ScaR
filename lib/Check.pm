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
			print "Step 0-1: Read length of fastq file has errors\n"; exit;
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
				($front, $back) = (split /\|/, $seq)[0, 1]; $front = substr($front, 1); # remove the first base of breakpoint sequence
			} elsif ( $seq =~/\*/ ) {
				($front, $back) = (split /\*/, $seq)[0, 1]; $back = substr($back, 0, length($back)-1); # remove the first base of breakpoint sequence
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
			$front_frag =~tr/[a-z]/[A-Z]/; $front_frag =~tr/ATCG/TAGC/; $front_frag=reverse($front_frag);
			$back_frag =~tr/[a-z]/[A-Z]/; $back_frag =~tr/ATCG/TAGC/; $back_frag=reverse($back_frag);
			$ref->{$header}[2] = [$front_tag, $front_frag]; # reverse complementary breakpoint seq of GeneA
			$ref->{$header}[3] = [$back_tag, $back_frag]; # reverse complementary breakpoint seq of GeneB
			#print "original: $ref->{$header}[0][1] $ref->{$header}[1][1]; reverse_complementary: $ref->{$header}[2][1] $ref->{$header}[3][1]\n"; #--testing
		}
		close IN; 
	}

	sub extract_coordinate {
		my ($gene, $geneA, $geneB, $read_l, $coord, $ref, $resource, $tran_seq, $genomic_seq, $hash_trans_seq_add) = @_;
		# // $gene: the path of Gene_hg38.txt
		# // $geneA: gene name of GeneA; // $geneB: gene name of GeneB
		# // $read_l: read length of paired-end fastq sequences
		# // $coord: genomic coordinates input format of --coordinate
		# // $ref: \%scaff_seq - main data structure
		# // $resource: "ensembl", "gencode" and "ucsc"
		# // $tran_seq: "ensembl_transcript.fa", "gencode_transcript.fa" and "ucsc_transcript.fa" - transcriptome sequences
		# // $genomic_seq: "GRCh38.primary_assembly.genome.fa" - genomic sequences
		# // $hash_trans_seq_add: collect genomic sequence as transcript references if breakpoint at intron or intergenic
		# separate each genomic breakpoint combination
		my @genomic_pos; # collect the given genomic breakpoint position
		foreach my $single (split /\,/, $coord) {
			$single =~s/^[\s]+//; $single =~s/[\s]+$//;
			my ($geneA_pos, $geneB_pos) = (split /\|/, $single)[0, 1]; 
			if ( defined($geneA_pos) and defined($geneB_pos) ) {
				my $chrom_geneA; if ( $geneA_pos =~/([\w]+)\:/ ) { $chrom_geneA = $1; $chrom_geneA =~s/[\s]+//; }
				my $chrom_geneB; if ( $geneB_pos =~/([\w]+)\:/ ) { $chrom_geneB = $1; $chrom_geneB =~s/[\s]+//; }
				my $coord_geneA; if ( $geneA_pos =~/\:([\d]+)/ ) { $coord_geneA = $1; $coord_geneA =~s/[\s]+//; }
				my $coord_geneB; if ( $geneB_pos =~/\:([\d]+)/ ) { $coord_geneB = $1; $coord_geneB =~s/[\s]+//; }
				my $strand_geneA; if ( $geneA_pos =~/\:(\+|\-)/ ) { $strand_geneA = $1; $strand_geneA =~s/[\s]+//; } else { $strand_geneA = undef; } # if strand not, undef
				my $strand_geneB; if ( $geneB_pos =~/\:(\+|\-)/ ) { $strand_geneB = $1; $strand_geneB =~s/[\s]+//; } else { $strand_geneB = undef; } # if strand not, undef
				push @genomic_pos, [$coord_geneA, $coord_geneB, $chrom_geneA, $chrom_geneB, $strand_geneA, $strand_geneB]; # breakpos of GeneA, breakpos of GeneB, Chrom of GeneA, Chrom of GeneB, Strand of GeneA, Strand of GeneB
			} else {
				print "\nStep 0-3: $single format for genomic coordinate of breakpoint is wrong, please set --coordinate in correct format\n\n"; exit;
			}
		}

		my %name; $/ = "\n"; # separator symbol
		my %chrom_include = ("1"=>1,"2"=>1,"3"=>1,"4"=>1,"5"=>1,"6"=>1,"7"=>1,"8"=>1,"9"=>1,"10"=>1,"11"=>1,"12"=>1,"13"=>1,"14"=>1,"15"=>1,"16"=>1,"17"=>1,"18"=>1,"19"=>1,"20"=>1,"21"=>1,"22"=>1,"X"=>1,"Y"=>1,"MT"=>1); # only consider the genes in chromosme (1..22, X, Y and MT); gene within scaffold is filtered out.
		open (IN, "cut -s -f1,5,10 $gene | sort | uniq |") || die "Step 0-3: Cannot open gene annotation file:$!\n";
		while ( <IN> ) {
			chomp $_; my ($ensembl, $symbol, $chrom) = (split /\t/, $_)[0, 1, 2]; 
			if ( exists($name{$symbol}) ) {
				if ( exists($chrom_include{$chrom}) ) {
					$name{$symbol} = $ensembl; # key => gene_symbol; value => ensembl_id
				}
			} else {
				$name{$symbol} = $ensembl; # key => gene_symbol; value => ensembl_id
			}
		}
		close IN;
		
		my %cdna; $/ = ">"; # separator symbol
		if ( $resource eq "ensembl" ) {
			open (IN, "$tran_seq") || die "Step 0-3: Cannot load ensembl transcript sequence dataset:$!\n";
			while ( <IN> ) {
				chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
				my @array = (split /\n/, $_); my $header = shift @array; $header =~s/\>//g;
				my $seq = join '', @array;
				my ($gene, $tran, $strand, $start, $end) = (split /\|/, $header)[0, 1, 2, 3, 4];
				my @first = split /\,/, $start; my @second = split /\,/, $end;
				push @{$cdna{$gene}}, [$seq, $strand, \@first, \@second, $tran]; # @first=(exon_start1, exon_start2, ...); @second=(exon_end1, exon_end2, ...);
			}
			close IN;
		} elsif ( $resource eq "gencode" ) {
			open (IN, "$tran_seq") || die "Step 0-3: Cannot load gencode transcript sequence dataset:$!\n";
			while ( <IN> ) {
				chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
				my @array = (split /\n/, $_); my $header = shift @array; $header =~s/\>//g;
				my $seq = join '', @array; 
				my ($gene, $tran, $strand, $start, $end) = (split /\|/, $header)[0, 1, 8, 9, 10];
				$tran =~s/\.[\d]+$//g; $gene =~s/\.[\d]+$//g;
				my @first = split /\,/, $start; my @second = split /\,/, $end;
				push @{$cdna{$gene}}, [$seq, $strand, \@first, \@second, $tran]; # @first=(exon_start1, exon_start2, ...); @second=(exon_end1, exon_end2, ...);
			}
			close IN;
		} else {
			my %name_tmp; # intermiddle file
			open (UCSC, "cut -s -f2,4,10,11,13 $resource |") || die "Step 0-3: Cannot load ucsc transcript annotation:$!\n";
			while ( <UCSC> ) {
				chomp $_; my ($tran, $strand, $start, $end, $gene) = (split /\t/, $_)[0, 1, 2, 3, 4]; $start =~s/\,$//; $end =~s/\,$//;
				my @first = split /\,/, $start; my @second = split /\,/, $end;
				if ( exists($name{$gene}) ) {
					$name_tmp{$tran} = [$name{$gene}, $strand, \@first, \@second, $tran]; # @first=(exon_start1, exon_start2, ...); @second=(exon_end1, exon_end2, ...);
				}
			}
			close UCSC;

			open (IN, "$tran_seq") || die "Step 0-3: Cannot load ucsc transcript sequence dataset:$!\n";
			while ( <IN> ) {
				chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
				my @array = (split /\n/, $_);
				my $header = shift @array; $header =~s/\>//g; my $tran = (split /\s/, $header)[0];
				if ( exists($name_tmp{$tran}) ) {
					my $seq = join '', @array;  $seq =~tr/atcg/ATCG/;
					push @{$cdna{$name_tmp{$tran}[0]}}, [$seq, $name_tmp{$tran}[1], $name_tmp{$tran}[2], $name_tmp{$tran}[3], $name_tmp{$tran}[4]]; # 
				}
			}
			close IN;
		}
 
		my %genomic_dna; # load genomic dna sequence
		open (IN, "$genomic_seq") || die "Step 0-3: Cannot load HG38 genomic sequence dataset:$!\n";
		while ( <IN> ) {
			chomp $_; next if (! defined($_) ); next if ( $_ eq '' );
			my @array = (split /\n/, $_);
			my $header = shift @array; $header =~s/\>//g; my $chrom_name = (split /\s/, $header)[0];
			my $seq = join '', @array; 
			$genomic_dna{$chrom_name} = $seq;
		}
		close IN;
		$/ = "\n";

		# judgement ensembl_id format or gene symbol format for --coordinate
		foreach my $ref_pos ( @genomic_pos ) { # // $ref_pos->[0]: breakpos of GeneA; # // $ref_pos->[1]: breakpos of GeneB
			my $x = 0; # index the number of combinations per each breakpoint pair
			my %scaff_seq_A; # collect scaffold breakpoint sequences of one given genomic coordinate for GeneA
			my %scaff_seq_B; # collect scaffold breakpoint sequences of one given genomic coordinate for GeneB
			my $strand_seq_A; # strand direction of GeneA
			my $strand_seq_B; # strand direction of GeneB
			if ( $geneA =~/ENSG/ ) { # geneA names as ensembl id
				if ( $geneB =~/ENSG/ )  { # geneB names as ensembl id
					if ( exists($cdna{$geneA}) ) { $strand_seq_A = &coordinate_seq($cdna{$geneA}, $read_l, $ref_pos->[0], "upstream", \%scaff_seq_A, $ref_pos->[4]); }
					if ( exists($cdna{$geneB}) ) { $strand_seq_B = &coordinate_seq($cdna{$geneB}, $read_l, $ref_pos->[1], "downstream", \%scaff_seq_B, $ref_pos->[5]); }
				} else { # geneB names as gene symbol
					if ( exists($cdna{$geneA}) ) { $strand_seq_A = &coordinate_seq($cdna{$geneA}, $read_l, $ref_pos->[0], "upstream", \%scaff_seq_A, $ref_pos->[4]); }
					if ( exists($name{$geneB}) ) {
						if ( exists($cdna{$name{$geneB}}) ) {
							$strand_seq_B = &coordinate_seq($cdna{$name{$geneB}}, $read_l, $ref_pos->[1], "downstream", \%scaff_seq_B, $ref_pos->[5]);
						} 
					} 
				}
			} else { # geneA names as gene symbol
				if ( $geneB =~/ENSG/ )  { # geneB names as ensembl id
					if ( exists($name{$geneA}) ) {
						if ( exists($cdna{$name{$geneA}}) ) {
							$strand_seq_A = &coordinate_seq($cdna{$name{$geneA}}, $read_l, $ref_pos->[0], "upstream", \%scaff_seq_A, $ref_pos->[4]);
						}
					}
					if ( exists($cdna{$geneB}) ) { $strand_seq_B = &coordinate_seq($cdna{$geneB}, $read_l, $ref_pos->[1], "downstream", \%scaff_seq_B, $ref_pos->[5]); }
				} else {
					if ( exists($name{$geneA}) ) {
						if ( exists($cdna{$name{$geneA}}) ) {
							$strand_seq_A = &coordinate_seq($cdna{$name{$geneA}}, $read_l, $ref_pos->[0], "upstream", \%scaff_seq_A, $ref_pos->[4]); 
						}
					}
					if ( exists($name{$geneB}) ) {
						if ( exists($cdna{$name{$geneB}}) ) {
							$strand_seq_B = &coordinate_seq($cdna{$name{$geneB}}, $read_l, $ref_pos->[1], "downstream", \%scaff_seq_B, $ref_pos->[5]);
						}
					}
				}
			}

			if ( %scaff_seq_A ) { # breakpoint site falls in exon region of GeneA, please use transcript sequence
				# remove redundency string in %scaff_seq_A
				my $size_A; my @long_A;
				foreach my $key ( reverse sort {length($a) <=> length($b)} keys %scaff_seq_A ) {
					if ( ! defined($size_A) ) {
						$size_A = length($key); push @long_A, $key;
					} else {
						if ( length($key) == $size_A ) { push @long_A, $key; }
					}
				}
				foreach my $key ( keys %scaff_seq_A ) {
					foreach my $sample ( @long_A ) {
						next if ( length($key) == length($sample) );
						if ( $sample =~/$key/ ) { delete($scaff_seq_A{$key}); }
					}
				}
			} else { # breakpoint site falls out of exon region of GeneA, please use genomic sequence
				if ( defined($strand_seq_A) ) {
					if ( exists($genomic_dna{$ref_pos->[2]}) ) { # if chromsome name matches for GeneA
						my $tmp_seq; # collect the breakpoint sequence
						my $tmp_genomic_seq; # collect the genomic sequence as transcript reference
						my $seed = 701; # set the length of genomic sequence as the transcript reference: 700-----|-----700
						my $tmp_seq_rev; # collect the breakpoint sequence (reverse complementary)
						my $tmp_genomic_seq_rev; # collect the genomic sequence (reverse complementary) as transcript reference
						my $seed_rev = 699; # set the length of genomic sequence (reverse complementary)as the transcript reference: 700-----|-----700
						# supporse geneA is at positive strand
						my $start_pos = $ref_pos->[0] - $read_l;
						$tmp_seq = substr($genomic_dna{$ref_pos->[2]}, $start_pos, $read_l);
						$start_pos = $ref_pos->[0] - $seed;
						$tmp_genomic_seq = substr($genomic_dna{$ref_pos->[2]}, $start_pos, $seed+$seed);
						# supporse geneA is at negative strand
						$tmp_seq_rev = substr($genomic_dna{$ref_pos->[2]}, $ref_pos->[0]-1, $read_l);
						$tmp_seq_rev =~tr/ATCG/TAGC/; $tmp_seq_rev = reverse($tmp_seq_rev); 
						$tmp_genomic_seq_rev = substr($genomic_dna{$ref_pos->[2]}, $ref_pos->[0]-1-$seed_rev, $seed_rev+$seed_rev);
						$tmp_genomic_seq_rev =~tr/ATCG/TAGC/; $tmp_genomic_seq_rev = reverse($tmp_genomic_seq_rev);

						if ( $strand_seq_A eq '+' ) { # if geneA at position strand
							if ( defined($ref_pos->[4]) ) { # if the strand direction input of coordinate for geneA present
								if ( $ref_pos->[4] eq $strand_seq_A ) {
									$scaff_seq_A{$tmp_seq} = [0, "noexon", $ref_pos->[0], "direction", 0-$seed, "trans$ref_pos->[0]pos", $tmp_genomic_seq];
								} else {
									$scaff_seq_A{$tmp_seq_rev} = [0, "noexon", $ref_pos->[0], "direction_rev", 0-$seed, "trans$ref_pos->[0]pos", $tmp_genomic_seq];
								}
							} else {
								$scaff_seq_A{$tmp_seq} = [0, "noexon", $ref_pos->[0], "direction", 0-$seed, "trans$ref_pos->[0]pos", $tmp_genomic_seq];
								$scaff_seq_A{$tmp_seq_rev} = [0, "noexon", $ref_pos->[0], "direction_rev", 0-$seed, "trans$ref_pos->[0]pos", $tmp_genomic_seq];
							}
						} else { # if geneA at negatvie strand
							if ( defined($ref_pos->[4]) ) { # if the strand direction input of coordinate for geneA present
								if ( $ref_pos->[4] eq $strand_seq_A ) {
									$scaff_seq_A{$tmp_seq_rev} = [0, "noexon", $ref_pos->[0], "direction", 0-$seed_rev, "trans$ref_pos->[0]neg", $tmp_genomic_seq_rev];
								} else {
									$scaff_seq_A{$tmp_seq} = [0, "noexon", $ref_pos->[0], "direction_rev", 0-$seed_rev, "trans$ref_pos->[0]neg", $tmp_genomic_seq_rev];
								}
							} else {
								$scaff_seq_A{$tmp_seq} = [0, "noexon", $ref_pos->[0], "direction_rev", 0-$seed_rev, "trans$ref_pos->[0]neg", $tmp_genomic_seq_rev];
								$scaff_seq_A{$tmp_seq_rev} = [0, "noexon", $ref_pos->[0], "direction", 0-$seed_rev, "trans$ref_pos->[0]neg", $tmp_genomic_seq_rev];
							}
						}
						$seed = $seed + 1;
						$seed_rev = $seed_rev - 1;
					} else {
						print "Step 0-3: Chromosome name incorrect for GeneA in --coordinate format\n"; exit;
					}
				} else {
					print "Step 0-3: $geneA is not in database, please input the valid gene\n"; exit;
				}
			}
			if ( %scaff_seq_B ) { # breakpoint site falls in exon region of GeneB, please use transcript sequence
				# remove redundency string in %scaff_seq_B
				my $size_B; my @long_B;
				foreach my $key ( reverse sort {length($a) <=> length($b)} keys %scaff_seq_B ) {
					if ( ! defined($size_B) ) {
						$size_B = length($key); push @long_B, $key;
					} else {
						if ( length($key) == $size_B ) { push @long_B, $key; }
					}
				}
				foreach my $key ( keys %scaff_seq_B ) {
					foreach my $sample ( @long_B ) {
						next if ( length($key) == length($sample) );
						if ( $sample =~/$key/ ) { delete($scaff_seq_B{$key}); }
					}
				}
			} else { # breakpoint site falls out of exon region of GeneB, please use genomic sequence
				if ( defined($strand_seq_B) ) {
					if ( exists($genomic_dna{$ref_pos->[3]}) ) { # if chromsome name matches for GeneB
						my $tmp_seq; # collect the breakpoint sequence
						my $tmp_genomic_seq; # collect the genomic sequence as transcript reference
						my $seed = 701; # set the length of genomic sequence as the transcript reference: 700-----|-----700
						my $tmp_seq_rev; # collect the breakpoint sequence (reverse complementary)
						my $tmp_genomic_seq_rev; # collect the genomic sequence (reverse complementary) as transcript reference
						my $seed_rev = 699; # set the length of genomic sequence (reverse complementary)as the transcript reference: 700-----|-----700
						# suppose that geneA at positive strand
						$tmp_seq = substr($genomic_dna{$ref_pos->[3]}, $ref_pos->[1]-1, $read_l);
						$tmp_genomic_seq = substr($genomic_dna{$ref_pos->[3]}, $ref_pos->[1]-1-$seed, $seed+$seed);
						# suppose that geneA at negative strand
						my $start_pos = $ref_pos->[1] - $read_l;
						$tmp_seq_rev = substr($genomic_dna{$ref_pos->[3]}, $start_pos, $read_l);
						$tmp_seq_rev =~tr/ATCG/TAGC/; $tmp_seq_rev = reverse($tmp_seq_rev);
						$start_pos = $ref_pos->[1] - $seed_rev;
						$tmp_genomic_seq_rev = substr($genomic_dna{$ref_pos->[3]}, $start_pos, $seed_rev+$seed_rev);
						$tmp_genomic_seq_rev =~tr/ATCG/TAGC/; $tmp_genomic_seq_rev = reverse($tmp_genomic_seq_rev);
	
						if ( $strand_seq_B eq '+' ) { # if geneB at position strand
							if ( defined($ref_pos->[5]) ) { # if the strand direction input of coordinate for geneB present
								if ( $ref_pos->[5] eq $strand_seq_B ) {
									$scaff_seq_B{$tmp_seq} = [0, "noexon", $ref_pos->[1], "direction", 0-$seed, "trans$ref_pos->[1]pos", $tmp_genomic_seq];
								} else {
									$scaff_seq_B{$tmp_seq_rev} = [0, "noexon", $ref_pos->[1], "direction_rev", 0-$seed, "trans$ref_pos->[1]pos", $tmp_genomic_seq];
								}
							} else {
								$scaff_seq_B{$tmp_seq} = [0, "noexon", $ref_pos->[1], "direction", 0-$seed, "trans$ref_pos->[1]pos", $tmp_genomic_seq];
								$scaff_seq_B{$tmp_seq_rev} = [0, "noexon", $ref_pos->[1], "direction_rev", 0-$seed, "trans$ref_pos->[1]pos", $tmp_genomic_seq];
							}
						} else {
							if ( defined($ref_pos->[5]) ) { # if the strand direction input of coordinate for geneB present
								if ( $ref_pos->[5] eq $strand_seq_B ) {
									$scaff_seq_B{$tmp_seq_rev} = [0, "noexon", $ref_pos->[1], "direction", 0-$seed_rev, "trans$ref_pos->[1]neg", $tmp_genomic_seq_rev];
								} else {
									$scaff_seq_B{$tmp_seq} = [0, "noexon", $ref_pos->[1], "direction_rev", 0-$seed_rev, "trans$ref_pos->[1]neg", $tmp_genomic_seq_rev];
								}
							} else {
								$scaff_seq_B{$tmp_seq} = [0, "noexon", $ref_pos->[1], "direction_rev", 0-$seed_rev, "trans$ref_pos->[1]neg", $tmp_genomic_seq_rev];
								$scaff_seq_B{$tmp_seq_rev} = [0, "noexon", $ref_pos->[1], "direction", 0-$seed_rev, "trans$ref_pos->[1]neg", $tmp_genomic_seq_rev];
							}
						}
						$seed = $seed + 1;
						$seed_rev = $seed_rev - 1;
					} else {
						print "Step 0-3: Chromosome name incorrect for GeneB in --coordinate format\n"; exit;
					}
				} else {
					print "Step 0-3: $geneB is not in database, please input the valid gene\n"; exit;
				}
			}

			# print all scaffold combination
			foreach my $seqA ( sort {$a cmp $b} keys %scaff_seq_A ) {
				foreach my $seqB ( sort {$a cmp $b} keys %scaff_seq_B ) {
					next if ( $scaff_seq_A{$seqA}[3] eq "direction_rev" and  $scaff_seq_B{$seqB}[3] eq "direction_rev" );
					my $header = "alt_".$scaff_seq_A{$seqA}[2]."_".$scaff_seq_B{$seqB}[2]."_".$x;
					$ref->{$header}[0] = [$scaff_seq_A{$seqA}[0], $seqA]; # breakpoint seq of GeneA
					$ref->{$header}[1] = [$scaff_seq_B{$seqB}[0], $seqB]; # breakpoint seq of GeneB
					
					my $rev_com_seqA = $seqA; $rev_com_seqA =~tr/[a-z]/[A-Z]/; $rev_com_seqA =~tr/ATCG/TAGC/; $rev_com_seqA=reverse($rev_com_seqA);
					my $rev_com_seqB = $seqB; $rev_com_seqB =~tr/[a-z]/[A-Z]/; $rev_com_seqB =~tr/ATCG/TAGC/; $rev_com_seqB=reverse($rev_com_seqB);
					$ref->{$header}[2] = [$scaff_seq_A{$seqA}[0], $rev_com_seqA]; # reverse complementary breakpoint seq of GeneA
					$ref->{$header}[3] = [$scaff_seq_B{$seqB}[0], $rev_com_seqB]; # reverse complementary breakpoint seq of GeneB
					# $ref->{$header}[4]/[5] = [gene_name, transcript_length, transcript_name, transcript_sequence]
					if ( $scaff_seq_A{$seqA}[1] eq "noexon" ) { $ref->{$header}[4] = [$geneA, $scaff_seq_A{$seqA}[4], $scaff_seq_A{$seqA}[5], $scaff_seq_A{$seqA}[6]]; }
					if ( $scaff_seq_B{$seqB}[1] eq "noexon" ) { $ref->{$header}[5] = [$geneB, $scaff_seq_B{$seqB}[4], $scaff_seq_B{$seqB}[5], $scaff_seq_B{$seqB}[6]]; }
					print "$header\t$scaff_seq_A{$seqA}[3].$scaff_seq_B{$seqB}[3]\t$scaff_seq_A{$seqA}[1]\t$scaff_seq_A{$seqA}[0]\t$seqA\t$scaff_seq_B{$seqB}[1]\t$scaff_seq_B{$seqB}[0]\t$seqB\n";
					$x++;
				}
			}
		}
	}

	sub coordinate_seq {
		my ($ref_seq, $read_range, $gene_pos, $type, $ref_hash, $given_strand) = @_;
		# // $ref_seq: transcript sequence array - $cdna{$name{$geneB}}
		# // $read_range: read length of paired-end fastq
		# // $gene_pos: breakpoint coordinate of GeneA / GeneB
		# // $type: "upstream" or "downstream"
		# // $ref_hash: $ref_hash->{"scaffold_breakpoint_sequence"} = [|Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_coordinate, "direction"/"direction_rev"]
		# // $given_strand: strand direction of geneA / geneB in coordinate input (i.e. undef, + or -)
		my $strand_direction; # set strand direction of gene partner
		foreach my $id ( @{$ref_seq} ) {
			$strand_direction = $id->[1];
			if ( $id->[1] eq '+' ) { # transcript in positive strand
				my $tag = 0; my $count = 0;
				for ( my $i=0; $i < scalar(@{$id->[2]}); $i++ ) {
					if ( $gene_pos >= $id->[2][$i] && $gene_pos <= $id->[3][$i] ) { # start < genomic_coordinate < end
						$count = $count + $gene_pos - $id->[2][$i] + 1;
						$tag = 1; last;
					} else {
						$count = $count + $id->[3][$i] - $id->[2][$i] + 1;
					}
				} 
				if ( $tag == 1 ) {
					if ( $type eq "upstream" ) { # upstream sequence
						my $front_frag; my $front_tag = 0;
						my $front_frag_reverse; my $front_tag_reverse = 0;
						if ( $count >= $read_range ) {
			 				$front_frag = substr($id->[0], $count-$read_range, $read_range);
						} else {
							$front_frag = substr($id->[0], 0, $count);
							$front_tag = $read_range - $count;
						}

						$front_frag_reverse = substr($id->[0], $count-1, $read_range); 
						$front_frag_reverse =~tr/ATCG/TAGC/; $front_frag_reverse = reverse($front_frag_reverse);
						my $count_back = length($id->[0]) - $count + 1;
						if ( $count_back < $read_range ) { $front_tag_reverse = $read_range - $count_back; }

						$front_frag = substr($front_frag, 1); $front_tag = $front_tag + 1; # get the length of upstream (read_length - 1)
						$front_frag_reverse = substr($front_frag_reverse, 1); $front_tag_reverse = $front_tag_reverse + 1; # get the length of upstream (read_length - 1)
						if ( defined($given_strand) ) {
							if ( $given_strand eq $id->[1] ) {
								$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							} else {
								$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							}
						} else {
							$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
						}
					} elsif ( $type eq "downstream" ) { # downstream sequence
						my $front_frag; my $front_tag = 0;
						my $front_frag_reverse; my $front_tag_reverse = 0;
						$front_frag = substr($id->[0], $count-1, $read_range);
						my $count_back = length($id->[0]) - $count + 1;
						if ( $count_back < $read_range ) { $front_tag = $read_range - $count_back; }

						if ( $count >= $read_range ) {
							$front_frag_reverse = substr($id->[0], $count-$read_range, $read_range);
						} else {
							$front_frag_reverse = substr($id->[0], 0, $count);
							$front_tag_reverse = $read_range - $count;
						}
						$front_frag_reverse =~tr/ATCG/TAGC/; $front_frag_reverse = reverse($front_frag_reverse);
						
						$front_frag = substr($front_frag, 0, length($front_frag)-1); $front_tag = $front_tag + 1; # get the length of upstream (read_length - 1)
						$front_frag_reverse = substr($front_frag_reverse, 0, length($front_frag_reverse)-1); $front_tag_reverse = $front_tag_reverse + 1; # get the length of upstream (read_length - 1)
						if ( defined($given_strand) ) {
							if ( $given_strand eq $id->[1] ) {
								$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_po
							} else {
								$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|,$transcript_id, genomic_pos
							}
						} else {
							$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
						}
					}
				}
			} else { # transcript in negative strand
				my $tag = 0; my $count = 0;
				for ( my $i=0; $i < scalar(@{$id->[2]}); $i++ ) {
					if ( $gene_pos >= $id->[2][$i] && $gene_pos <= $id->[3][$i] ) { # start < genomic_coordinate < end
						$count = $count + $id->[3][$i] - $gene_pos + 1;
						$tag = 1; last;
					} else {
						$count = $count + $id->[3][$i] - $id->[2][$i] + 1;
					}
				}
				if ( $tag == 1 ) {
					if ( $type eq "upstream" ) { # upstream sequence
						my $front_frag; my $front_tag = 0;
						my $front_frag_reverse; my $front_tag_reverse = 0;

						if ( $count >= $read_range ) {
							$front_frag = substr($id->[0], $count-$read_range, $read_range);
						} else {
							$front_frag = substr($id->[0], 0, $count);
							$front_tag = $read_range - $count;
						}

						$front_frag_reverse = substr($id->[0], $count-1, $read_range);
						$front_frag_reverse =~tr/ATCG/TAGC/; $front_frag_reverse = reverse($front_frag_reverse);
						my $count_back = length($id->[0]) - $count + 1;
						if ( $count_back < $read_range ) { $front_tag_reverse = $read_range - $count_back; }

						$front_frag = substr($front_frag, 1); $front_tag = $front_tag + 1; # get the length of upstream (read_length - 1)
						$front_frag_reverse = substr($front_frag_reverse, 1); $front_tag_reverse = $front_tag_reverse + 1; # get the length of upstream (read_length - 1)
						if ( defined($given_strand) ) {
							if ( $given_strand eq $id->[1] ) {
								$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							} else {
								$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							}
						} else {
							$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
						}
					} elsif ( $type eq "downstream" ) { # downstream sequence
						my $front_frag; my $front_tag = 0;
						my $front_frag_reverse; my $front_tag_reverse = 0;

						$front_frag = substr($id->[0], $count-1, $read_range);
						my $count_back = length($id->[0]) - $count + 1;
						if ( $count_back < $read_range ) { $front_tag = $read_range - $count_back; }

						if ( $count >= $read_range ) {
							$front_frag_reverse = substr($id->[0], $count-$read_range, $read_range);
						} else {
							$front_frag_reverse = substr($id->[0], 0, $count);
							$front_tag_reverse = $read_range - $count;
						}
						$front_frag_reverse =~tr/ATCG/TAGC/; $front_frag_reverse = reverse($front_frag_reverse);

						$front_frag = substr($front_frag, 0, length($front_frag)-1); $front_tag = $front_tag + 1; # get the length of upstream (read_length - 1)
						$front_frag_reverse = substr($front_frag_reverse, 0, length($front_frag_reverse)-1); $front_tag_reverse = $front_tag_reverse + 1; # get the length of upstream (read_length - 1)
						if ( defined($given_strand) ) {
							if ( $given_strand eq $id->[1] ) {
								$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							} else {
								$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							}
						} else {
							$ref_hash->{$front_frag} = [$front_tag, $id->[4], $gene_pos, "direction"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
							$ref_hash->{$front_frag_reverse} = [$front_tag_reverse, $id->[4], $gene_pos, "direction_rev"]; # |Exp(read_length)-Obs(read_length)|, $transcript_id, genomic_pos
						}
					}
				}
			}
		}
		return($strand_direction);
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
			if ( exists($name{$symbol}) ) {
				if ( exists($chrom_include{$chrom}) ) {
					$name{$symbol} = $ensembl; # key => gene_symbol; value => ensembl_id
				}
			} else {
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
				my ($gene, $tran) = (split /\|/, $header)[0, 1];
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
				my ($gene, $tran) = (split /\|/, $header)[0, 1]; $tran =~s/_//g;
				if (! defined($gene) ) { print "Header of user-defined transcript references are error, please follow >gene|transcript\n"; exit; }
				if (! defined($tran) ) { print "Header of user-defined transcript references are error, please follow >gene|transcript\n"; exit; }
				my $seq = join '', @array; my $len = length($seq);
				$cdna{$gene}{$len} = [$tran, $seq];
			}
			close IN;
		} else {
			print "Step1: User-defind transcript reference sequence not present, please continue\n";
			if ( defined($scaffold->[4]) ) { # if genomic sequences as transcript references for GeneA as when using --coordinate
				$cdna{$scaffold->[4][0]}{$scaffold->[4][1]} = [$scaffold->[4][2], $scaffold->[4][3]];
			}
			if ( defined($scaffold->[5]) ) { # if genomic sequences as transcript references for GeneB as when using --coordinate
				$cdna{$scaffold->[5][0]}{$scaffold->[5][1]} = [$scaffold->[5][2], $scaffold->[5][3]];
			}
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
					print "Step 1: $geneA is valid gene\n";	my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
					$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
				} else {
					print "Step 1: $geneA is not in database, please input the valid gene or use user-defined sequences\n";
				}

				if ( exists($cdna{$geneB}) ) {
					print "Step 1: $geneB is valid gene\n";	my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
					$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
				} else {
					print "Step 1: $geneB is not in database, please input the valid gene or use user-defined sequences\n";
				}
			} else { # geneB names as gene symbol
				if (exists($cdna{$geneA}) ) {
					print "Step 1: $geneA is valid gene\n"; my @array; $scaff_final{$geneA} = \@array; # define the sub- data structure
					$tag_scaff_A = &output_geneA($cdna{$geneA}, $scaffold, $scaff_final{$geneA});
				} else {
					print "Step 1: $geneA is not in database, please input the valid gene or use user-defined sequences\n";
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
					print "Step 1: $geneB is valid gene\n"; my @array; $scaff_final{$geneB} = \@array; # define the sub- data structure
					$tag_scaff_B = &output_geneB($cdna{$geneB}, $scaffold, $scaff_final{$geneB});
				} else {
					print "Step 1: $geneB is not in database, please input the valid gene or use user-defined sequences\n";
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
			my @cdna_geneB = split /$split_tag/, $cdna_seq->{$rank}[1];  #print ">$cdna_seq->{$rank}[0]\n$cdna_seq->{$rank}[1]\n";
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
				@cdna_geneB = split /$split_tag/, $cdna_seq->{$rank}[1]; #print "$split_tag\n", scalar(@cdna_geneB), "\n";
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
					push @{$scaff_seq->[0]}, $ref->[0]; push @{$scaff_seq->[1]}, $ref->[1]; push @{$scaff_seq->[2]}, $ref->[2]; # [sequence, transcript_id, break_seq_geneB]
				}
			}
			return(1);
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
					push @{$scaff_seq->[0]}, $ref->[0]; push @{$scaff_seq->[1]}, $ref->[1]; push @{$scaff_seq->[2]}, $ref->[2]; # [sequence, transcript_id, break_seq_geneA]
				}
			}
			return(1);
		} else {
			return(0);
		}
	}
1;
