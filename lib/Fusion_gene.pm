package Fusion_gene;

use strict;
use warnings;

	sub cDNA {
		my ($ref, $path, $a, $b) = @_;
		my $tag_a = 0; my $tag_b = 0;
		my %ensembl; my %gene_name;
		my %chrom_include = ("1"=>1,"2"=>1,"3"=>1,"4"=>1,"5"=>1,"6"=>1,"7"=>1,"8"=>1,"9"=>1,"10"=>1,"11"=>1,"12"=>1,"13"=>1,"14"=>1,"15"=>1,"16"=>1,"17"=>1,"18"=>1,"19"=>1,"20"=>1,"21"=>1,"22"=>1,"X"=>1,"Y"=>1,"MT"=>1); # only consider the genes in chromosme (1..22, X, Y and MT); gene within scaffold is filtered out.
		open (IN, "cut -s -f1,3,4,5,9,10 $path | sort | uniq |") || die "Step 3: cannot open this path for cDNA path:$!\n";
		while ( <IN> ) {
			chomp $_; my ($gene_en, $start, $end, $symbol, $type, $chr) = (split /\t/, $_)[0,1,2,3,4,5];
			if (! defined($gene_en) ) { $gene_en = ""; }
			if (! defined($symbol) ) { $symbol = ""; }

			if ( exists($gene_name{$symbol}) ) {
				if ( exists($chrom_include{$chr}) ) {
					$gene_name{$symbol} = [$chr, $start, $end];
				}
			} else {
				$gene_name{$symbol} = [$chr, $start, $end];
			}
			$ensembl{$gene_en} = [$chr, $start, $end];
		}
		close IN;

		if ( exists($ensembl{$a}) ) {
			push @{$ref}, [$ensembl{$a}[0], $ensembl{$a}[1], $ensembl{$a}[2], $a]; $tag_a = 1;
		} elsif ( exists($gene_name{$a}) ) {
			push @{$ref}, [$gene_name{$a}[0], $gene_name{$a}[1], $gene_name{$a}[2], $a]; $tag_a = 1;
		}

		if ( exists($ensembl{$b}) ) {
			push @{$ref}, [$ensembl{$b}[0], $ensembl{$b}[1], $ensembl{$b}[2], $b]; $tag_b = 1;
		} elsif ( exists($gene_name{$b}) ) {
			push @{$ref}, [$gene_name{$b}[0], $gene_name{$b}[1], $gene_name{$b}[2], $b]; $tag_b = 1;
		}

		return($tag_a, $tag_b);
	}
1;
