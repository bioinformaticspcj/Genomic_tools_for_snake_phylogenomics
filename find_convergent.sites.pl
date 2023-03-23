#!/usr/bin/perl -w 
use strict;
if($#ARGV< 3){
	print " 
	pcj finished 2021.1.14 
	This script can tackle 3 kinds of variations( all->one,one->one,one->deletion);
	Implement provean with NR vertebrate database 2020.7 
	usage:
	perl $0 Gshe,Dacu,Pbiv aligned.fas -2.5(score threshold) Smer(represented background_species) 1> out 2> log\n";
	exit;
}
my @targets=split/,/,$ARGV[0];
my $file=$ARGV[1];
my $back_ground=$ARGV[3];
my $back_ground_candinate;
my %hash_target;
my %hash_target1;
my %hash_others;
my %others;
my %other_seq;
$/=">";
my $score_threshold=$ARGV[2] ||=-2.5;
for(@targets){
	$hash_target1{$_}=();
}
open IN,"$file" or die "can not open $file\n";
while(<IN>){
	chomp;
	next unless($_);
	my @lines=split/\n/;
	my @names=split/[_\|]/,$lines[0];
	my $seq=join"",@lines[1..$#lines];
	$back_ground_candinate=$names[0];
	$other_seq{$names[0]}=$seq if(!exists $hash_target1{$names[0]});
	my @aa=split"",$seq;
	push @{$hash_target{$names[0]}},@aa if(exists $hash_target1{$names[0]});
	push @{$hash_others{$names[0]}},@aa if(!exists $hash_target1{$names[0]});
}


my %tar_sites;
my $target_exact;
my @gaps;
my %tar_gap=();
foreach my $tar (@targets){
	if( exists $hash_target{$tar} ){
		$target_exact++;
		for(my $i=0;$i<=$#{$hash_target{$tar}};$i++){
			push @{$tar_sites{$i}},${$hash_target{$tar}}[$i];
			$gaps[$i]=0 if(${$hash_target{$tar}}[$i] ne "-");
			$tar_gap{$i}++,$gaps[$i]=${$hash_target{$tar}}[$i] if(${$hash_target{$tar}}[$i] eq "-");
		}
	}
}
my $file_out=0;
my $out;
foreach my $po(keys %tar_sites){
	#print "$po====\n";
	my %count;
	my @tar_get_site=grep{++$count{$_}<2} @{$tar_sites{$po}};
#	print @{$tar_sites{"1866"}},"========\n";
 	if($#tar_get_site==0 && $tar_get_site[0] ne "-"){
		my @other_sites;
		my $gap=0;
		my $flag=0;
		foreach my $other(keys %hash_others){
			if(${$hash_others{$other}}[$po] eq $tar_get_site[0]){
				$flag=1;
				next;
			}elsif(${$hash_others{$other}}[$po] ne $tar_get_site[0] ){
				push @other_sites,${$hash_others{$other}}[$po];	
				$gap++ if(${$hash_others{$other}}[$po] eq "-");
			}		
		}
		if($flag==0){
			my %count;
			my $po1=$po+1;
			my @all_other_sites=grep{++$count{$_}<2} @other_sites;
			if($#all_other_sites==0 || ($#all_other_sites==1 && $gap>=1 && $gap/($#other_sites+1) <=0.5)){
				$file_out=1;
				$out.="$po1:$all_other_sites[0]->$tar_get_site[0]:";
				$out.="strict\t";
			}elsif($#all_other_sites> 0){
				$file_out=1;
				$out.="$po1:unknow->$tar_get_site[0]:";
				$out.="loose\t";
			}
		}
	}
	if($#tar_get_site==0 && $tar_get_site[0] eq "-"){
		my @other_sites1;
		my @other_sites2;
		my $gap=0;
		my $flag=0;
		foreach my $other(keys %hash_others){
                        if(${$hash_others{$other}}[$po] ne $tar_get_site[0] ){
                                push @other_sites1,${$hash_others{$other}}[$po]; 
                        }else{        
				push @other_sites2,${$hash_others{$other}}[$po];
                                $gap++ if(${$hash_others{$other}}[$po] eq "-");
			}
		}
		my %count; 
                my $po1=$po+1;
                my @all_other_sites1=grep{++$count{$_}<2} @other_sites1;
		if($#all_other_sites1==0 && $gap>=1 && $gap/($#other_sites1+$#other_sites2+2) <=0.5){
			$file_out=1;
			$out.="$po1:$all_other_sites1[0]->$tar_get_site[0]:";
			$out.="strict999\t";	
		}
	}	
	if($#tar_get_site==1 && $gaps[$po] eq "-" && $tar_gap{$po}/$target_exact <0.3){
		my $rate=$tar_gap{$po}/$target_exact;
		my $j;
		for(my$i=0;$i<=$#tar_get_site;$i++){
			$j=$i if ($tar_get_site[$i] eq "-");
		}	
		my @other_sites;
                my $gap=0;
                my $flag=0;
		my $po1=$po+1;
                foreach my $other(keys %hash_others){
                	if(${$hash_others{$other}}[$po] eq $tar_get_site[1-$j]){
                        	$flag=1;
                                last;
                        }else{
                        	push @other_sites,${$hash_others{$other}}[$po];
                                $gap=1 if(${$hash_others{$other}}[$po] eq "-");
                        }
		}
                if($flag==0){
                	my %count;
                        my @all_other_sites=grep{++$count{$_}<2} @other_sites;
			my $c=1;
                	for(my$i=0;$i<=$#all_other_sites;$i++){
                        	$c=$i if ($all_other_sites[$i] eq "-");
                	}

                        if($#all_other_sites==0 || ($#all_other_sites==1 && $gap==1)){
                        	$file_out=1;
                                $out.="$po1:$all_other_sites[1-$c]->$tar_get_site[1-$j]:" if($#all_other_sites==1);
				$out.="$po1:$all_other_sites[0]->$tar_get_site[1-$j]:" if($#all_other_sites==0 && $all_other_sites[0] ne "-");
                                $out.="strict:$rate\t" if($#all_other_sites==1 || ($#all_other_sites==0 && $all_other_sites[0] ne "-"));
                        }elsif($#all_other_sites> 0){
                        	$file_out=1;
                                $out.="$po1:unknow->$tar_get_site[1-$j]:";
                                $out.="loose\t";
                        }
		}
		
	}
					
}
 if($out){
	my @sites=split/\t/,$out;
	my @aa;
	my @old_position;
	@aa=split"",$other_seq{$back_ground} if(exists $other_seq{$back_ground});
	@aa=split"",$other_seq{$back_ground_candinate} if(!exists $other_seq{$back_ground});
	print STDERR"no $back_ground found ! select $back_ground_candinate as back ground !\n" if(!exists $other_seq{$back_ground});
	
	open OUT, ">$file.var" or die "can not open $file.var\n";
	foreach my $site(@sites){
		my @infor=split/:/,$site;
		my ($old,$new)=split/\->/,$infor[1];
		if($infor[2] eq "strict"){
			$aa[$infor[0]-1]=$old if($aa[$infor[0]-1] eq "-" && $old ne "-");
		}
	}
	foreach my $site(@sites){
		my $gaps=0;
		my @infor=split/:/,$site;
		foreach my $aa_site(@aa[0..$infor[0]-2]){
			$gaps++ if($aa_site eq "-");
		}
                my ($old,$new)=split/\->/,$infor[1]; 
		my $po=$infor[0]-$gaps;
		print OUT "$old$po$new\n" if($old ne "-" && $infor[2] eq "strict" && $new ne "-");
		print OUT "$old${po}del\n" if($old ne "-" && $infor[2] eq "strict" && $new eq "-");	
		push @old_position,$infor[0] if($old ne "-" && $infor[2] eq "strict" && $new ne "-");
		push @old_position,$infor[0] if($old ne "-" && $infor[2] eq "strict" && $new eq "-");
	}
	close OUT;
	my $other_seq1=join"",@aa;
	$other_seq1=~s/\-//g;
	open OUT,">$file.other.seq.fas" or die "can not open $file.other.seq.fas\n";
	print OUT">background\n$other_seq1\n";
	close OUT;
#=head
	`rm -rf $file.var.tmp log` if(-e "$file.var.tmp");
	`mkdir $file.var.tmp`;
	#print "/data/nfs/OriginTools/varitions_predicts/provean-1.1.5_build/bin/provean.sh -q $file.other.seq.fas -v $file.var --save_supporting_set $file.other.seq.fas.sss --tmp_dir tmp --num_threads 50 2>log |grep -v \"#\"\n";
#=head
	my $result=`/data/nfs/OriginTools/varitions_predicts/provean-1.1.5_build/bin/provean.sh -q $file.other.seq.fas -v $file.var --save_supporting_set $file.other.seq.fas.sss --tmp_dir $file.var.tmp --num_threads 100 2>>log |grep -v "#"|grep -v "\\["`;
	#print "/data/nfs/OriginTools/varitions_predicts/provean-1.1.5_build/bin/provean.sh -q $file.other.seq.fas -v $file.var --save_supporting_set $file.other.seq.fas.sss --tmp_dir $file.var.tmp --num_threads 100 2>>log |grep -v \n";
	my @results=split/\n/,$result;
	for(my $i=0;$i<=$#results;$i++){
		chomp $results[$i];
		#print "$re===$score_threshold\n";
		my ($site,$score)=split" ",$results[$i];
		if($score=~/[A-Za-z]/){
			print STDERR "error: /data/nfs/OriginTools/varitions_predicts/provean-1.1.5_build/bin/provean.sh -q $file.other.seq.fas -v $file.var --save_supporting_set $file.other.seq.fas.sss --tmp_dir $file.var.tmp --num_threads 100 2>>log no site found!\n";
			exit;
		}
		$out.="$old_position[$i]:$site:$score:pass\t" if($score <= $score_threshold);			
	}
#=cut
	$out=~s/\t$//;	
}
unlink "$file.other.seq.fas";
unlink "$file.var";
`rm -rf $file.var.tmp`;

print "$file\t$out\n" if($file_out==1 && $out);
