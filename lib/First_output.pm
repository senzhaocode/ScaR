package First_output;
use strict;
use warnings;

	sub grep_spanning {
		my ($dis_ref, $fastq_1, $fastq_2, $path, $sep) = @_;
		open (OUT, ">>$path/read_mapped_info") || die "Step 2-3: cannot open summary output path:$!\n";
		open (OUT1, ">$path/spanning_1.txt") || die "Step 2-3: cannot open first round spanning_1 output path:$!\n";
		open (OUT2, ">$path/spanning_2.txt") || die "Step 2-3: cannot open second round spanning_2 output path:$!\n";
		foreach my $id ( keys %{$dis_ref} ) {
			my $header = "@".$id.$sep;

			if ( scalar(@{$dis_ref->{$id}}) == 2 ) {
				my $seq_1_tr = $dis_ref->{$id}[0][3]; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
				my $seq_1_qual = $dis_ref->{$id}[0][6]; $seq_1_qual = reverse($seq_1_qual);
				my $seq_2_tr = $dis_ref->{$id}[1][3]; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
				my $seq_2_qual = $dis_ref->{$id}[1][6]; $seq_2_qual = reverse($seq_2_qual);
				#/* use flag value for filtering
				if ( $dis_ref->{$id}[0][5] == 81 ) { # judge the flag of first element = 1st end (reverse complement)
					#if ( $dis_ref->{$id}[1][5] == 161 ) { # judge the flag of second element = 2nd end
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2]\n"; 
						print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
					#} 
				} elsif ( $dis_ref->{$id}[0][5] == 97 ) { # judge the flag of first element = 1st end
					#if ( $dis_ref->{$id}[1][5] == 145 ) { # judge the flag of second element = 2nd end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2](reverse complement)\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
					#}
				} elsif ( $dis_ref->{$id}[0][5] == 65 ) { # judge the flag of first element = 1st end
					#if ( $dis_ref->{$id}[1][5] == 129 ) { # judge the flag of second element = 2nd end
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2]\n"; 
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n"; 
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
					#}
				} elsif ( $dis_ref->{$id}[0][5] == 113 ) { # judge the flag of first element = 1st end (reverse complement)
					#if ( $dis_ref->{$id}[1][5] == 177 ) { # judge the flag of second element = 2nd end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2](reverse_complement)\n";
						print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n"; 
						print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
					#}
				} elsif ( $dis_ref->{$id}[0][5] == 161 ) { # judge the flag of first element = 2nd end
					#if ( $dis_ref->{$id}[1][5] == 81 ) { # judge the flag of second element = 1st end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2]\n";
						print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
					#}
				} elsif ( $dis_ref->{$id}[0][5] == 145 ) { # judge the flag of first element = 2nd end (reverse complement)
					#if ( $dis_ref->{$id}[1][5] == 97 ) { # judge the flag of second element = 1st end
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](reverse complement)\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
						print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
					#}
				} elsif ( $dis_ref->{$id}[0][5] == 129 ) { # judge the flag of first element = 2nd end
					#if ( $dis_ref->{$id}[1][5] == 65 ) { # judge the flag of second element = 1st end
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2]\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
					#}
				} elsif ( $dis_ref->{$id}[0][5] == 177 ) { # judge the flag of first element = 2nd end (reverse complement)
					#if ( $dis_ref->{$id}[1][5] == 113 ) { # judge the flag of second element = 1st end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2](reverse_complement)\n";
						print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
						print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
					#}
				} else {
					if ( $dis_ref->{$id}[0][5] == 83 ) { # judge the flag of first element = 1st end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2]\n";
						print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
					} elsif ( $dis_ref->{$id}[0][5] == 99 ) { # judge the flag of first element = 1st end
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2](reverse complement)\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
					} elsif ( $dis_ref->{$id}[0][5] == 163 ) { # judge the flag of first element = 2nd end
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2]\n";
						print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
					} elsif ( $dis_ref->{$id}[0][5] == 147 ) { # judge the flag of first element = 2nd end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](reverse complement)\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
						print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
					} elsif ( $dis_ref->{$id}[0][5] == 67 ) { # judge the flag of first element = 1st end
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2]\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
					} elsif ( $dis_ref->{$id}[0][5] == 115 ) { # judge the flag of first element = 1st end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2](reverse_complement)\n";
						print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
						print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
					} elsif ( $dis_ref->{$id}[0][5] == 131 ) { # judge the flag of first element = 2nd end
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2]\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
					} elsif ( $dis_ref->{$id}[0][5] == 179 ) { # judge the flag of first element = 2nd end (reverse complement)
						print OUT "spanning\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2](reverse_complement)\n";
						print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
						print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
					} else {
						print "spanning\t$id has incorrect flag value\n";
					}
				}
			}
		}
		close OUT1;
		close OUT2;
		close OUT;
	}

	sub grep_discordant {
		my ($dis_ref, $fastq_1, $fastq_2, $path, $sep) = @_;
		open (OUT, ">$path/read_mapped_info") || die "Step 2-3: cannot open summary output path:$!\n";
        	open (OUT1, ">$path/discordant_split_1.txt") || die "Step 2-3: cannot open first round discordant_1 output path:$!\n";
        	open (OUT2, ">$path/discordant_split_2.txt") || die "Step 2-3: cannot open second round discordant_2 output path:$!\n";

		foreach my $id ( keys %{$dis_ref} ) {
			my $header = "@".$id.$sep;

			if ( scalar(@{$dis_ref->{$id}}) == 2 ) {
				my $seq_1_tr = $dis_ref->{$id}[0][3]; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
				my $seq_1_qual = $dis_ref->{$id}[0][6]; $seq_1_qual = reverse($seq_1_qual);
				my $seq_2_tr = $dis_ref->{$id}[1][3]; $seq_2_tr =~tr/ATCG/TAGC/; $seq_2_tr = reverse($seq_2_tr);
				my $seq_2_qual = $dis_ref->{$id}[1][6]; $seq_2_qual = reverse($seq_2_qual);

				if ( $dis_ref->{$id}[0][2] eq "scaffold" && $dis_ref->{$id}[1][2] eq "scaffold") { # both two end mapped to scaffold
					if ( $dis_ref->{$id}[0][5] == 99 ) { # judge the flag of first element = 1st end
						if ( $dis_ref->{$id}[1][5] == 147 ) { # judge the flag of second element = 2nd end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2](reverse complement)\n";
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
							print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 83 ) { # judge the flag of first element = 1st end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 163 ) { # judge the flag of second element = 2nd end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse complement)\t$dis_ref->{$id}[1][2]\n";
							print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 67 ) { # judge the flag of first element = 1st end
						if ( $dis_ref->{$id}[1][5] == 131 ) { # judge the flag of second element = 2nd end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2]\n";
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 115 ) { # judge the flag of first element = 1st end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 179 ) { # judge the flag of second element = 2nd end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse complement)\t$dis_ref->{$id}[1][2](reverse complement)\n";
							print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
							print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 147 ) { # judge the flag of first element = 2nd end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 99 ) { # judge the flag of second element = 1st end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](reverse complement)\n";
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
							print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 163 ) { # judge the flag of first element = 2nd end
						if ( $dis_ref->{$id}[1][5] == 83 ) { # judge the flag of second element = 1st end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](reverse complement)\t$dis_ref->{$id}[0][2]\n";
							print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 131 ) { # judge the flag of first element = 2nd end
						if ( $dis_ref->{$id}[1][5] == 67 ) { # judge the flag of second element = 1st end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2]\n";
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 179 ) { # judge the flag of first element = 2nd end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 115 ) { # judge the flag of second element = 1st end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](reverse complement)\t$dis_ref->{$id}[0][2](reverse complement)\n";
							print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
							print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
						}
					} else {
						print "discordant_split\t$id has incorrect flag value\n";
					}
				} else {
					if ( $dis_ref->{$id}[0][5] == 81 ) { # judge the flag of first element = 1st end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 161 ) { # judge the flag of second element = 2nd end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2]\n";
							print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 97 ) { # judge the flag of first element = 1st end
						if ( $dis_ref->{$id}[1][5] == 145 ) { # judge the flag of second element = 2nd end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2](reverse complement)\n";
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
							print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 65 ) { # judge the flag of first element = 1st end
						if ( $dis_ref->{$id}[1][5] == 129 ) { # judge the flag of second element = 2nd end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[1][2]\n"; 
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
 						}
					} elsif ( $dis_ref->{$id}[0][5] == 113 ) { # judge the flag of first element = 1st end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 177 ) { # judge the flag of second element = 2nd end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[1][2](reverse_complement)\n";
							print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
							print OUT2 "@"."$id/2\n$seq_2_tr\n+\n$seq_2_qual\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 161 ) { # judge the flag of first element = 2nd end
						if ( $dis_ref->{$id}[1][5] == 81 ) { # judge the flag of second element = 1st end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2]\n";
							print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 145 ) { # judge the flag of first element = 2nd end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 97 ) { # judge the flag of second element = 1st end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2](reverse complement)\n";
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
							print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 129 ) { # judge the flag of first element = 2nd end
						if ( $dis_ref->{$id}[1][5] == 65 ) { # judge the flag of second element = 1st end
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2]\t$dis_ref->{$id}[0][2]\n";
							print OUT1 "@"."$id/1\n$dis_ref->{$id}[1][3]\n+\n$dis_ref->{$id}[1][6]\n";
							print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						}
					} elsif ( $dis_ref->{$id}[0][5] == 177 ) { # judge the flag of first element = 2st end (reverse complement)
						if ( $dis_ref->{$id}[1][5] == 113 ) { # judge the flag of second element = 1nd end (reverse complement)
							print OUT "discordant_split\t$id\t$dis_ref->{$id}[1][2](reverse_complement)\t$dis_ref->{$id}[0][2](reverse_complement)\n";
							print OUT1 "@"."$id/1\n$seq_2_tr\n+\n$seq_2_qual\n";
							print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
						}
					} else {
						print "discordant_split\t$id has incorrect flag value\n";
					}
                        	}
			} elsif ( scalar(@{$dis_ref->{$id}}) == 1 ) { # Not included 
				if ( $dis_ref->{$id}[0][2] eq "scaffold" ) {
					my $seq_1_tr = $dis_ref->{$id}[0][3]; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
					my $seq_1_qual = $dis_ref->{$id}[0][6]; $seq_1_qual = reverse($seq_1_qual);
					my $hit_1; my $hit_2; # set judgement

					if ( $dis_ref->{$id}[0][5] == 99 ) {
						print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2]\t$dis_ref->{$id}[0][4](reverse_complement)\n";
						print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
						if ( $fastq_2 =~/\.gz/ ) { # compressed
							$hit_2 = `zgrep -A 3 "$header" $fastq_2`; chomp $hit_2;
						} else { # uncompressed
							$hit_2 = `grep "$header" -A 3 $fastq_2`; chomp $hit_2;
						}
						my @array = (split /\n/, $hit_2);
						print OUT2 "@"."$id/2\n$array[1]\n+\n$array[3]\n";
					} elsif ( $dis_ref->{$id}[0][5] == 147 ) {
						print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][4]\t$dis_ref->{$id}[0][2](reverse_complement)\n";
						if ( $fastq_1 =~/\.gz/ ) { # compressed
							$hit_1 = `zgrep -A 3 "$header" $fastq_1`; chomp $hit_1;
						} else { # uncompressed
							$hit_1 = `grep "$header" -A 3 $fastq_1`; chomp $hit_1;
						}
						my @array = (split /\n/, $hit_1);
						print OUT1 "@"."$id/1\n$array[1]\n+\n$array[3]\n";
						print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
					} elsif ( $dis_ref->{$id}[0][5] == 83 ) {
						print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\t$dis_ref->{$id}[0][4]\n";
						print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
						if ( $fastq_2 =~/\.gz/ ) { # compressed
							$hit_2 = `zgrep -A 3 "$header" $fastq_2`; chomp $hit_2;
						} else { # uncompressed
							$hit_2 = `grep "$header" -A 3 $fastq_2`; chomp $hit_2;
						}
						my @array = (split /\n/, $hit_2);
						print OUT2 "@"."$id/2\n$array[1]\n+\n$array[3]\n";
					} elsif ( $dis_ref->{$id}[0][5] == 163 ) {
						print OUT "discordant_split\t$id\t$dis_ref->{$id}[0][4](reverse_complement)\t$dis_ref->{$id}[0][2]\n";
						if ( $fastq_1 =~/\.gz/ ) { # compressed
							$hit_1 = `zgrep -A 3 "$header" $fastq_1`; chomp $hit_1;
						} else { # uncompressed
							$hit_1 = `grep "$header" -A 3 $fastq_1`; chomp $hit_1;
						}
						my @array = (split /\n/, $hit_1);
						print OUT1 "@"."$id/1\n$array[1]\n+\n$array[3]\n";
						print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][6]\n";
					} else {
						print "discordant_split\t$id has incorrect flag value\n";
					}
				}
			}
		}
		close OUT1;
		close OUT2;
		close OUT;
	}
	
	sub grep_singlton {
		my ($dis_ref, $fastq_1, $fastq_2, $path, $sep) = @_;
		open (OUT, ">>$path/read_mapped_info") || die "Step 2-3: cannot open summary output path:$!\n";
		open (OUT1, ">$path/singlton_split_1.txt") || die "Step 2-3: cannot open first round singlton_1 output path:$!\n";
		open (OUT2, ">$path/singlton_split_2.txt") || die "Step 2-3: cannot open first round singlton_2 output path:$!\n";
	
		foreach my $id ( keys %{$dis_ref} ) {
			my $header = "@".$id.$sep;

			if ( scalar(@{$dis_ref->{$id}}) == 1 ) {
				my $hit_1; my $hit_2; # set judgement
				my $seq_1_tr = $dis_ref->{$id}[0][3]; $seq_1_tr =~tr/ATCG/TAGC/; $seq_1_tr = reverse($seq_1_tr);
				my $seq_1_qual = $dis_ref->{$id}[0][5]; $seq_1_qual = reverse($seq_1_qual);

				if ( $dis_ref->{$id}[0][4] == 73 ) {
					print OUT "singlton_split\t$id\t$dis_ref->{$id}[0][2]\tNULL\n"; 
					print OUT1 "@"."$id/1\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][5]\n";
					if ( $fastq_2 =~/\.gz/ ) { # compressed
						$hit_2 = `zgrep -A 3 "$header" $fastq_2`; chomp $hit_2;
					} else {
						$hit_2 = `grep "$header" -A 3 $fastq_2`; chomp $hit_2;
					}
					my @array = (split /\n/, $hit_2);
					print OUT2 "@"."$id/2\n$array[1]\n+\n$array[3]\n";
				} elsif ( $dis_ref->{$id}[0][4] == 89 || $dis_ref->{$id}[0][4] == 121 ) {
					print OUT "singlton_split\t$id\t$dis_ref->{$id}[0][2](reverse_complement)\tNULL\n"; 
					print OUT1 "@"."$id/1\n$seq_1_tr\n+\n$seq_1_qual\n";
					if ( $fastq_2 =~/\.gz/ ) { # compressed
						$hit_2 = `zgrep -A 3 "$header" $fastq_2`; chomp $hit_2;
					} else {
						$hit_2 = `grep "$header" -A 3 $fastq_2`; chomp $hit_2;
					}
					my @array = (split /\n/, $hit_2);
					print OUT2 "@"."$id/2\n$array[1]\n+\n$array[3]\n";
				} elsif ( $dis_ref->{$id}[0][4] == 137 ) {
					print OUT "singlton_split\t$id\tNULL\t$dis_ref->{$id}[0][2]\n"; 
					if ( $fastq_1 =~/\.gz/ ) { # compressed
						$hit_1 = `zgrep -A 3 "$header" $fastq_1`; chomp $hit_1;
					} else {
						$hit_1 = `grep "$header" -A 3 $fastq_1`; chomp $hit_1;
					}
					my @array = (split /\n/, $hit_1);
					print OUT1 "@"."$id/1\n$array[1]\n+\n$array[3]\n";
					print OUT2 "@"."$id/2\n$dis_ref->{$id}[0][3]\n+\n$dis_ref->{$id}[0][5]\n";
				} elsif ( $dis_ref->{$id}[0][4] == 153 || $dis_ref->{$id}[0][4] == 185 ) {
					print OUT "singlton_split\t$id\tNULL\t$dis_ref->{$id}[0][2](reverse_complement)\n"; 
					if ( $fastq_1 =~/\.gz/ ) { # compressed
						$hit_1 = `zgrep -A 3 "$header" $fastq_1`; chomp $hit_1;
					} else {
						$hit_1 = `grep "$header" -A 3 $fastq_1`; chomp $hit_1;
					}
					my @array = (split /\n/, $hit_1);
					print OUT1 "@"."$id/1\n$array[1]\n+\n$array[3]\n";
					print OUT2 "@"."$id/2\n$seq_1_tr\n+\n$seq_1_qual\n";
				} else {
					print "singlton_split\t$id has incorrect flag value\n";
				}
			}
		}
		close OUT1;
		close OUT2;
		close OUT;
	}
1;
