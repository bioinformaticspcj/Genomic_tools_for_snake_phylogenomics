#!/usr/bin/perl
use strict;
use warnings;
if($#ARGV < 0){
	print "$0 Bed.fa > merged.CNE.fa\n";
	exit; 
}
my %seqs;
opendir DIR,"$ARGV[0]" or die "can not open $ARGV[0]\n";
while(readdir DIR){
	next,unless(/\.final/);
	my @CNE_fa_files=glob"$ARGV[0]/$_/*.fa",;
	foreach my $file(@CNE_fa_files){
		next,unless(-e "$file.anc.dnd");
		open FA,"$file" or die "can not open $file\n";
		$/=">";
		#print STDERR "open $file \n";
		while(<FA>){
			chomp;	
			next,unless($_);
			my @lines = split/\n/;
			my @names = split/\./,$lines[0];
		
			$seqs{$names[0]}.=join"",@lines[1..$#lines];
		}
		$/="\n";	
	}
}
foreach(keys %seqs){
	my $out_seq=&format_seq($seqs{$_});
	print">$_\n$out_seq\n";
}
sub format_seq{
	my $seq = shift;
	my $length = length $seq;
	my $left_num = $length % 50; 
	my @seqs = $seq =~/.{50}/g;
	my $out = join"\n",@seqs;
	my $left_seq = substr($seq,$length-$left_num,$left_num);
	$out = $out."\n$left_seq";
	return "$out";
}
