#!/usr/bin/perl
use FindBin qw($Bin $RealBin);
die "perl $0 <pep>\n" if $ARGV<0;
my $pep=shift;
my $pep_prank="$pep.best.fas";
my $pep_prank_trimal="$pep.best.fas.trimal";
my $pep_prank_html="$pep.best.fas.trimal.html";

#my $pep_prank_2="$pep.pep.prank.fasta.mark.2";
#my $pep_2="$pep.2";
    #}
`/data/tools/prank-msa/bin/prank -d=$pep -o=$pep`;
`/data/tools/trimal/source/trimal  -in  $pep_prank  -automated1 -out $pep_prank_trimal  -fasta  -htmlout  $pep_prank_html  `;
