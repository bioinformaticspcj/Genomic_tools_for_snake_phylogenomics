#!/usr/bin/perl -w
use strict;
#use Parallel::ForkManager;
my $forground=$ARGV[0];
my $background=$ARGV[1];
#my $max_process =$ARGV[3] ||4;
#my $pm = new Parallel::ForkManager($max_process);
my $all=1;
$all=$ARGV[2] if(defined $ARGV[2]);
my $p_adjust_methods="BH" ||shift;
if(!$forground || !$background){
	print " 
taget.txt(targeted gene names|IDs,one name|ID per line) 
background(g:profiler GMT format file) 
1(optional 1:for annotated genes only default:1) 0:for all input genes;
2(optional select second column to enrich default: 1);
example: 
perl $0 taget.txt test.GMT 1\n";
	exit;		
}
my $sep=$ARGV[3] ||1;
open TARGET,"$forground" or die "can not open $forground\n";
my (@background,%forground,@p_value,@p_ajust,@ano_genes);
while(<TARGET>){
	chomp;
	my @lines=split/[\t,]/;
	$forground{uc($lines[$sep-1])}=1;
}
close TARGET;

my $f_total=keys %forground;
my $ano_total_f=0;
my $ano_total_b=0;
my %ano_total_f;
my %ano_total_b;
open BACK,"$background" or die "can not open $background\n";
while(<BACK>){
	chomp;
	my @lines=split/\t/;
	my $ano_f=0;
	my $ano_b=$#lines-2+1;
	my @for_genes;
	my $flag=0;
	foreach my $gene_in_back (@lines[2..$#lines]){
		$ano_total_b{$gene_in_back}=1;
		if(exists $forground{uc($gene_in_back)}){
			$ano_total_f{$gene_in_back}=1;
			$ano_f++;
			$flag=1;	
			push @for_genes,$gene_in_back;	
		}
	}
	if($flag){
		my $back_gene=join",",@lines[2..$#lines];
		my $for_gene=join",",@for_genes;
		my $infor=join"===",(@lines[0..1],$ano_f,$ano_b,$for_gene,$back_gene);
		push @background,$infor;
	}
}
close  BACK;
if($#background<0){
	print "no forground gene found in background !\n";
	exit;
}
my $ano_total_f1=keys %ano_total_f;
my $ano_total_b1=keys %ano_total_b;
#LINE:
foreach my $item (@background){
#	$pm->start and next LINE;
	my @infor=split/===/,$item;
	my $pvalue;
	$pvalue=&get_over_presented($infor[2],$ano_total_f1,$infor[3],$ano_total_b1) if($all);
	$pvalue=&get_over_presented($infor[2],$f_total,$infor[3],$ano_total_b1) if(!$all);
	push @p_value,$pvalue;
#	$pm->finish;
	
}
#$pm->wait_all_children;

my $pvalue=join",",@p_value;

my $R_script=<<"RSCRIPT";
	p.ajust<-p.adjust(c($pvalue),"BH");
	print(p.ajust);
RSCRIPT
open OUT,">p.ajust.R" or die "can not open p.ajust.R\n";
print OUT "$R_script";
close OUT;
my $p_ajust_raw=`Rscript p.ajust.R`;
unlink "p.ajust.R";
chomp $p_ajust_raw;
@p_ajust=split" ",$p_ajust_raw;
@p_ajust=grep{!/\[/} @p_ajust;
print"term ID\tdescription\ttarget in path\tback_ground in path\ttarget genes\tback_ground genes\ttarget total\tback_ground total\tP-value\tP.adjust\n";
for(my $i=0;$i<=$#background;$i++){
	$background[$i]=~s/===/\t/g;
	print "$background[$i]\t$ano_total_f1\t$ano_total_b1\t$p_value[$i]\t$p_ajust[$i]\n" if($all);
	print "$background[$i]\t$f_total\t$ano_total_b1\t$p_value[$i]\t$p_ajust[$i]\n" if(!$all);
}
sub get_over_presented{
	my ($for_count,$for_total_count,$back_count,$back_total_count)=@_;
	my $for=$for_count-1;
	my $R_script=<<"RSCRIPT";
	pvalue<-phyper($for,$for_total_count,$back_total_count-$for_total_count,$back_count,lower.tail=F);
	print(pvalue)
RSCRIPT
	my $pvalue=`echo "$R_script" | Rscript -`;
	chomp $pvalue;
	$pvalue=~s/\[1\]\s+//;
	return $pvalue;
}
