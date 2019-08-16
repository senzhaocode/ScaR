package Simulate;
use strict;
use warnings;

	sub load {
		my ($simulate_ref, $simulate_ref_ensembl, $synthetic_path) = @_;
		open (IN, "$synthetic_path") || die "can not open synthetic fusion transcripts path:$!\n";
		while ( <IN> ) {
			chomp $_; my $line = $_; next if ( $line =~/Ensembl_id/ );
			my ($geneA, $enA, $chromA, $posA, $strandA, $geneB, $enB, $chromB, $posB, $strandB) = (split /\t/, $line)[0, 1, 2, 3, 5, 7, 8, 9, 10, 12];
			$simulate_ref->{$geneA}{$geneB} = [$strandA, $posA, $strandB, $posB];
			$simulate_ref_ensembl->{$enA}{$enB} = [$strandA, $posA, $strandB, $posB];
		}
		close IN;
	}
1;
