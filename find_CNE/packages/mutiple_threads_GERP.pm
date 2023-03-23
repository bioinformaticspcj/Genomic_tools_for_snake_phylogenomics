#!/usr/bin/perl
package mutiple_threads_GERP; 
require Exporter; 
my @ISA= qw(Exporter); 
my @EXPORT  = qw(block_file_deal); 
my @version = 1.0; 
use strict;
use warnings;
use Parallel::ForkManager;
sub block_file_deal{
	my($maf,$tree,$id,$seqname,$ref_species,$block,$threads,$cup)=@_;
	my @block_files=glob"$maf.block*.mfa";
	my $max_process= int ($cup/$threads);
	print "Threads for GERP: $max_process\n";
	my $pm = new Parallel::ForkManager($max_process);
	LINE:
	foreach my $block_file (@block_files){
		$pm->start and next LINE;
		if(-e $block_file){
			print "/data/nfs/OriginTools/call_CNE_new/GERP/gerpcol -t $tree -f $block_file -x .GERP.rates -e $ref_species ===\n";
                	system("/data/nfs/OriginTools/call_CNE_new/GERP/gerpcol -t $tree -f $block_file -x .GERP.rates -e $ref_species");
                	unlink "$maf";
                	system("/data/nfs/OriginTools/call_CNE_new/GERP/gerpelem -w .ex -f $block_file.GERP.rates -d 0.01");
                	open ELE,"$block_file.GERP.rates.elems" or die "can not open $block_file.GERP.rates.elems \n";
			open Tele,">$block_file.GERP.rates.elems.f" or die "can not open $block_file.GERP.rates.elems.f \n";
                	my $block_num=$1,if($block_file=~/(block\d+)/);
                	while(<ELE>){
                        	chomp;
                        	my $left=0;
                        	my $right=0;
                        	my @line = split/\t+/;
                        	$line[2]=$line[2]+$block->{$block_num};
                        	$line[1]=$line[1]+$block->{$block_num};
                        	my $out1 = "$seqname"."\t".$line[1]."\t".$line[2]."\t"."$id.GERP\t0\t+";
                        	print Tele "$out1\n",if($out1);
                	}
                	close ELE;
			close Tele;
		}
		else{exit;}
		$pm->finish;	
	}
	$pm->wait_all_children;
	system("cat $maf*.block*.GERP.rates.elems.f > ../CNE/$maf.c.con.bed.cname.GERP.bed");
	system("perl -p -i -ne 's/GERP/GERP\$./' ../CNE/$maf.c.con.bed.cname.GERP.bed");
}
1;
__END__ 
