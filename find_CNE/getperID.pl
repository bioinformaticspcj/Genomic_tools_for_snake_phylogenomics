#!/usr/bin/perl 
use strict;
use warnings;
use File::Spec;
my $path_curf = File::Spec->rel2abs(__FILE__);
my ($vol, $dir, $file) = File::Spec->splitpath($path_curf);
##### run prank #####
`$ARGV[1] -d=$ARGV[0] -o=$ARGV[0] -showtree -showanc -keep -prunetree -seed=10  `;
print "$ARGV[1] -d=$ARGV[0] -o=$ARGV[0] -showtree -showanc -keep -prunetree -seed=10\n";
#### parse tree ####
`Rscript $dir/parse_tree.R $ARGV[0].anc.dnd`;
`cat $ARGV[0].anc.dnd.* > $ARGV[0].alltree_infor`;
#### hash fasta #####
open FA,"$ARGV[0].anc.fas" or die "can not open $ARGV[0].anc.fas";
my %hash_fa;
$/=">";
while(<FA>){
	chomp;
	next,unless($_);
	my @line=split/\n/;
	my $head=$1,if($line[0]=~/(^\S+)/);
	$head=~s/^(\S+?)_.*/$1/g;
	my $seq=join"",@line[1..$#line];
	$hash_fa{$head}=$seq;
}
$/="\n";
close FA;
##### struct tree hash ####
my %parents_edge=();
my %edges_tag=();
my @children=();
my @parents=();
open TREE_INFOR,"$ARGV[0].alltree_infor" or die "can not open $ARGV[0].alltree_infor\n";
while(<TREE_INFOR>){
	chomp;
	my @line=split/\t/;
	if($_=~/^\d/){
		$parents_edge{$line[0]}.=$line[1].'-';
	}
	if($_=~/^[A-Z]/){
		$line[0]=~s/(\S+?)_.*/$1/g;
		push @children,$line[0];
	#	print "$_++\n";	
	}
	if($_=~/^#/){
		push @parents,$line[0];
		
	}
}
close TREE_INFOR;
@children=sort @children;
my $number_child=scalar @children;
my $i=0;
#=cut
my %parents_tag=();
my ($id1,$id2);
foreach my $p(keys %parents_edge){
	my($left,$right)=split/\-/,$parents_edge{$p};
	my $son1 = $left > $number_child ? $parents[$left-1-$number_child] : $children[$left-1-$number_child];
	my $son2 = $right > $number_child ? $parents[$right-1-$number_child] : $children[$right-1-$number_child];
	$parents_tag{$parents[$p-1-$number_child]} = "$son1-$son2";
	#print "$son1-$son2\n";
}
%parents_tag=&tag_parents(%parents_tag);
sub tag_parents{
	my %pa =@_;
	foreach (keys %pa){
		if($pa{$_}=~/(#\d+#)\-/ || $pa{$_}=~/\-(#\d+#)/){
			my $tag=$1;
			$pa{$_}=~s/$tag/$pa{$tag}/;
			#print "$pa{$_}\n";
			#&tag_parents(%pa),if($pa{$_}=~/#/);
		}
	}
	my @values= values %pa;
	my $str=join"",@values;
	%pa=&tag_parents(%pa),if($str=~/#/);
	return %pa;
}
#foreach (keys %parents_tag){
#	print "$_:$parents_tag{$_}\n";
#}

#print "$number_child====\n";
#my $a=keys %parents_edge;
#print "$a\n";
#print "$_\n" for @parents;
#=cut
open PERIDLOCAL,">$ARGV[0].perlocalid";
open ALLPERID,">$ARGV[0].perglobleid";
open ALLPERID_H,"> perglobleid_head"; 
my @cneid=split/\//,$ARGV[0];
$cneid[-1]=~s/\.fa$//;
#print "$cneid[-1]====\n";
my $id = $cneid[-1];
#my$id=$id[-1];

foreach ( sort{$a<=>$b} keys %parents_edge){
	$i++;	
	#print $parents[$i-1],"\n";
	my($left,$right)=split/\-/,$parents_edge{$_};
	if($left > $number_child){
		$id1 = &getID("$parents[$left-$number_child-1]",$hash_fa{$parents[$i-1]},$hash_fa{$parents[$left-$number_child-1]},$id);
	}else{
	 	$id1 = &getID("$children[$left-1]",$hash_fa{$parents[$i-1]},$hash_fa{$children[$left-1]},$id);
		#$parents_tag{$parents[$_-$number_child-1]}.= $children[$left-1]."-";
	 
	}
	if($right > $number_child){
	 	$id2 =  &getID("$parents[$right-$number_child-1]",$hash_fa{$parents[$i-1]},$hash_fa{$parents[$right-$number_child-1]},$id) ;
	}else{
	 	$id2 = &getID("$children[$right-1]",$hash_fa{$parents[$i-1]},$hash_fa{$children[$right-1]},$id);
		#$parents_tag{$parents[$_-$number_child-1]}.= $children[$right-1]."-";
	}
	#print "$parents[$_-$number_child-1]\t$parents_tag{$parents[$_-$number_child-1]}===\n",if(exists $parents_tag{$parents[$_-$number_child-1]});
	$id1=~s/$parents[$left-$number_child-1]/$parents_tag{$parents[$left-$number_child-1]}/g,if($parents[$left-$number_child-1]);
	$id2=~s/$parents[$right-$number_child-1]/$parents_tag{$parents[$right-$number_child-1]}/g,if($parents[$right-$number_child-1]);
	#print "$id1\t$id2\n$parents[$_-$number_child-1]\t$parents_tag{$parents[$_-$number_child-1]}\n";
	print PERIDLOCAL "$id1\n$id2\n";
	#print "$id1\n$id2\n";
}
my @species;
my @ids;
foreach my $child(@children){
	my $id = &getID("$child",$hash_fa{$parents[0]},$hash_fa{$child},$id);
#	print "$child\t$hash_fa{$parents[0]}\t$hash_fa{$child}\t$id===\n";
	my @out=split/\t/,$id;
	push @species,$child;
	push @ids,$out[2];
	#print "$child\n";
}
my $out_species=join"\t",@species;
print ALLPERID_H "species\t$out_species\n";
my $out_ids=join"\t",@ids;
print ALLPERID "$id\t$out_ids\n";
close PERIDLOCAL;
close ALLPERID;
close ALLPERID_H;

sub getID{
	my ($species,$ref,$qury,$id)=@_;
	my @seq1=split"",$ref;
	my @seq2=split"",$qury;
	my $len=length($ref);
	my $same=0;
	for(my $i=0;$i<=$#seq1;$i++){
		if($seq1[$i] eq $seq2[$i] && $seq1[$i] =~/[ATCGatcg]/ && $seq2[$i] =~/[ATCGatcg]/){
			$same++;
		}
	}
	#print "$species\t$same\t$len\t====$ref\t$qury\n";
	my $identity;
	eval { $identity=$same/$len};
	if($@){
		print STDERR "$ref\n$qury\n something wrong !\n";
		exit;
	}else{	
		return "$species\t$id\t$identity\t$len";
	}
}
