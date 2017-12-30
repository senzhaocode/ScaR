package Fusion_gene;

use strict;
use warnings;

	sub cDNA {
		my ($ref, $path, $a, $b) = @_;
		my $tag_a = 0; my $tag_b = 0;
		my %chrom_include = ("1"=>1,"2"=>1,"3"=>1,"4"=>1,"5"=>1,"6"=>1,"7"=>1,"8"=>1,"9"=>1,"10"=>1,"11"=>1,"12"=>1,"13"=>1,"14"=>1,"15"=>1,"16"=>1,"17"=>1,"18"=>1,"19"=>1,"20"=>1,"21"=>1,"22"=>1,"X"=>1,"Y"=>1,"MT"=>1); # only consider the genes in chromosme (1..22, X, Y and MT); gene within scaffold is filtered out.
		open (IN, "cut -s -f1,3,4,5,9,10 $path | sort | uniq |") || die "Step 3: cannot open this path for cDNA path:$!\n";
		while ( <IN> ) {
			chomp $_; my ($gene_en, $start, $end, $symbol, $type, $chr) = (split /\t/, $_)[0,1,2,3,4,5];
			next if (! exists($chrom_include{$chr}) );

			if ( $a eq $gene_en ) {
				push @{$ref}, [$chr, $start, $end, $a]; $tag_a = 1;
			} elsif ( $a eq $symbol ) {
				push @{$ref}, [$chr, $start, $end, $a]; $tag_a = 1;
			}

			if ( $b eq $gene_en ) {
				push @{$ref}, [$chr, $start, $end, $b]; $tag_b = 1;
			} elsif ( $b eq $symbol ) {
				push @{$ref}, [$chr, $start, $end, $b]; $tag_b = 1;
			}
		}
		close IN;
		if ( $tag_a == 1 && $tag_b == 1 ) {
		} else {
			print "Step 3: Gene name does not match to annotation\n"; exit;
		}
	}
1;
