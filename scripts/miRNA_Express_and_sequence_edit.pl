#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2011/5/20
#Modified: Gong Jing 2011-12-6
#Description: solexa miRNA express and sequence
my $version=1.00;

use strict;
use Getopt::Long;

my ($infile)=@ARGV;
my (%opts, $outdir);
open (OPT,"<$infile") or die "cannot open config.properties\n";


while(my $line=<OPT>){
                if($line=~/RESULT_DIRECTORY\s*=\s*(\S+)/)     { $opts{'out'}=$1;}
		        if($line=~/GENE_LIST_FILE\s*=\s*(\S+)/)       { $opts{'g'}=$1;}
		        if($line=~/SPECIES_ALIAS_NAME\s*=\s*(\S+)/)   { $opts{'species'}=$1;}
				}
				close OPT;
				$opts{'i'}=$opts{'out'}."\/prediction\/predictions";
				$opts{'p'}=$opts{'out'}."\/prediction\/precursors.fa"; 
			    $opts{'list'}=$opts{'out'}."\/target/microRNA.list";
				$opts{'fa'}=$opts{'out'}."\/target/microRNA.fa";

my $targetfile=$opts{'out'}."\/target/";
if( ! -e $targetfile ){ 
`mkdir "$targetfile"` ;
}
my $tag=$opts{'species'};
my $filein=$opts{'i'};
my $fileout=$opts{'list'};
my $out=$opts{'fa'};

my %hash_pri;
open PRI,"<$opts{p}";
while (my $aline=<PRI>) {
	chomp $aline;
	if($aline=~/^>(\S+)/){$hash_pri{$1}=$aline;}
}
close PRI;
my %gene;
if (defined $opts{'g'}) {
	open G,"<$opts{'g'}";
	while (my $aline=<G>) {
		chomp $aline;
		next if($aline=~/^\#/);
		my @temp=split/\t/,$aline;
		$temp[1]=~/\.(\S+)/;
		push @{$gene{$1}},$temp[2]."\t".$temp[3]."\t".$temp[4];
	}
	close G;
}

open IN,"<$filein"; #input file  
open OUT,">$fileout"; #output file  
open FA ,">$out";
print OUT "#ID_1\ttag\tchrom\tstrand\tpre_start\tpre_end\tm_strat\tm_end\ts_start\ts_end\tmature\tstar\tpre\tid_2\tm_exp\ts_exp\tm_exp_total\ts_exp_total\n";
my ($novel,$conserve,%uniq_id);
while (my $aline=<IN>) {
	chomp $aline;
	my $signal="Conserve";
	until ($aline =~ /^score_nucleus/){
		$aline = <IN>;
		if (eof) {last;}
	}
	if (eof) {last;}
#%%%%%%%%%%%%%%%%%%%%%%%% conserve %%%%%%%%%%%%%%%%%%%%%%5	
	if ($aline=~/^score_nucleus	3/) {
		$conserve++;
########## miRNA ID ################
		$aline=<IN>;
		chomp $aline;
		my @temp=split/\t/,$aline;
		my $mature=join ";",@temp;
		my @tmp=split/-/,$temp[0];
		$tmp[0]=$tag;
		my $id=join "-",@tmp;#my $first_id=$id;
		if (!defined $uniq_id{$id}) {
	$uniq_id{$id}=0;
	print OUT "$id\t";
	print FA ">$id\n";
		}
		else{
	$uniq_id{$id}++;
	print OUT "$id","_$uniq_id{$id}\t";
	print FA ">$id","_$uniq_id{$id}\n";
		}
########### annotate####################
		do {$aline=<IN>;} until($aline=~/flank_first_end/) ;
		chomp $aline;
		my @flank1=split/\t/,$aline;
		do {$aline=<IN>;} until($aline=~/flank_second_beg/) ;
		chomp $aline;
		my @flank2=split/\t/,$aline;
#		
########## mature start loop pre ####
		do {$aline=<IN>;} until($aline=~/mature_beg/) ;
		chomp $aline;
		my @start=split/\t/,$aline;
		$start[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/mature_end/) ;
		chomp $aline;
		my @end=split/\t/,$aline;
		$end[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/mature_seq/) ;
		chomp $aline;
		my @arr1=split/\t/,$aline;
		do {$aline=<IN>;} until($aline=~/pre_seq/) ;
		chomp $aline;
		my @arr2=split/\t/,$aline;
		do {$aline=<IN>;} until($aline=~/pri_id/) ;
		chomp $aline;
		my @pri_id=split/\t/,$aline;
my $annotation=&annotate($flank1[1],$flank2[1],$pri_id[1]);
		do {$aline=<IN>;} until($aline=~/star_beg/) ;
		chomp $aline;
		my @star_start=split/\t/,$aline;
		$star_start[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/star_end/) ;
		chomp $aline;
		my @star_end=split/\t/,$aline;
		$star_end[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/star_seq/) ;
		chomp $aline;
		my @arr3=split/\t/,$aline;
#		my $substring=substr($arr1[1],1,8);
#		$signal="Known" if($hash_known{$first_id} && ($hash_known{$first_id} =~/$substring/));

		print OUT "$signal\t$annotation\t$start[1]\t$end[1]\t$star_start[1]\t$star_end[1]\t";
		print OUT "$arr1[1]\t$arr3[1]\t$arr2[1]\t$mature\t";
		print FA "$arr1[1]\n";
########## reads count #############
		<IN>;
		my $count1=0;my $count2=0;my $count3=0;my $count4=0;
		$aline=<IN>;
		do {
	chomp $aline;
	my @reads=split/\t/,$aline;
	my @pos=();
	$reads[5]=~/(\d+)\.\.(\d+)/;
	$pos[0] =$1-$flank1[1];
	$pos[1] =$2-$flank1[1];
	$reads[0]=~/_x(\d+)$/;
#	$count3 +=$1 if($pos[0]<=$end[1] && $pos[1]>=$start[1] );
#	$count4 +=$1 if($pos[0]<=$star_end[1] && $pos[1]>=$star_start[1] );
#	$count1 =$1 if($pos[0]<=$end[1] && $pos[1]>=$start[1] && $count1<$1);
#	$count2 =$1 if($pos[0]<=$star_end[1] && $pos[1]>=$star_start[1] && $count2<$1);
	$count3 +=$1 if($end[1]-$pos[0]>=10 && $pos[1]-$start[1]>=10 );
	$count4 +=$1 if($star_end[1]-$pos[0]>=10 && $pos[1]-$star_start[1]>=10 );
	$count1 =$1 if($end[1]-$pos[0]>=10 && $pos[1]-$start[1]>=10 && $count1<$1);
	$count2 =$1 if($star_end[1]-$pos[0]>=10 && $pos[1]-$star_start[1]>=10 && $count2<$1);
	$aline=<IN>;
	chomp $aline;
		} until(length $aline < 1) ;
		print OUT "$count1\t$count2\t$count3\t$count4\n";
	}
#%%%%%%%%%%%%%%%%%%%%%%%% novel %%%%%%%%%%%%%%%%%%%%%%5	
	else{
		$novel++;
		$signal="Novel";
########### annotate####################
		do {$aline=<IN>;} until($aline=~/flank_first_end/) ;
		chomp $aline;
		my @flank1=split/\t/,$aline;
		do {$aline=<IN>;} until($aline=~/flank_second_beg/) ;
		chomp $aline;
		my @flank2=split/\t/,$aline;
#		
########## mature start loop pre ####
		do {$aline=<IN>;} until($aline=~/mature_beg/) ;
		chomp $aline;
		my @start=split/\t/,$aline;
		$start[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/mature_end/) ;
		chomp $aline;
		my @end=split/\t/,$aline;
		$end[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/mature_seq/) ;
		chomp $aline;
		my @arr1=split/\t/,$aline;
		do {$aline=<IN>;} until($aline=~/pre_seq/) ;
		chomp $aline;
		my @arr2=split/\t/,$aline;
		do {$aline=<IN>;} until($aline=~/pri_id/) ;
		chomp $aline;
		my @pri_id=split/\t/,$aline;
my $annotation=&annotate($flank1[1],$flank2[1],$pri_id[1]);
		do {$aline=<IN>;} until($aline=~/star_beg/) ;
		chomp $aline;
		my @star_start=split/\t/,$aline;
		$star_start[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/star_end/) ;
		chomp $aline;
		my @star_end=split/\t/,$aline;
		$star_end[1] -=$flank1[1];
		do {$aline=<IN>;} until($aline=~/star_seq/) ;
		chomp $aline;
		my @arr3=split/\t/,$aline;
		print OUT "$tag-miR-c-$novel\t$signal\t$annotation\t$start[1]\t$end[1]\t$star_start[1]\t$star_end[1]\t";
		print OUT "$arr1[1]\t$arr3[1]\t$arr2[1]\t\/\t";
		print FA ">$tag-miR-c-$novel\n$arr1[1]\n";
########## reads count #############
		<IN>;
		my $count1=0;my $count2=0;my $count3=0;my $count4=0;
		$aline=<IN>;
		do {
	chomp $aline;
	my @reads=split/\t/,$aline;
	my @pos=();
	$reads[5]=~/(\d+)\.\.(\d+)/;
	$pos[0] =$1-$flank1[1];
	$pos[1] =$2-$flank1[1];
	$reads[0]=~/_x(\d+)$/;
#	$count3 +=$1 if($pos[0]<=$end[1] && $pos[1]>=$start[1] );
#	$count4 +=$1 if($pos[0]<=$star_end[1] && $pos[1]>=$star_start[1]);
#	$count1 =$1 if($pos[0]<=$end[1] && $pos[1]>=$start[1] && $count1<$1);
#	$count2 =$1 if($pos[0]<=$star_end[1] && $pos[1]>=$star_start[1] && $count2<$1);
	$count3 +=$1 if($end[1]-$pos[0]>=10 && $pos[1]-$start[1]>=10 );
	$count4 +=$1 if($star_end[1]-$pos[0]>=10 && $pos[1]-$star_start[1]>=10 );
	$count1 =$1 if($end[1]-$pos[0]>=10 && $pos[1]-$start[1]>=10 && $count1<$1);
	$count2 =$1 if($star_end[1]-$pos[0]>=10 && $pos[1]-$star_start[1]>=10 && $count2<$1);
	$aline=<IN>;
	chomp $aline;
		} until(length $aline < 1) ;
		print OUT "$count1\t$count2\t$count3\t$count4\n";
	}
}

close IN;
close OUT;

sub annotate{
	my ($st,$ed,$id)=@_;
	my @temp=split/\s+/,$hash_pri{$id};
	my $ret;
	$id=~/_\d+/;
	my $chrom=$`;
	$ret .= $chrom."\t";
	$temp[1]=~/:/;
	my $strand=$';
	$ret .=$strand."\t";
	$temp[2]=~/:(\d+)/;
	my $start=$1+$st;
	$ret.= $start . "\t";
#	$temp[3]=~/:(\d+)/;
	my $end=$start+$ed-$st-2;
	$ret .=$end;
	if($gene{$chrom}){
		foreach my $gene (@{$gene{$chrom}}) {
	my @tmp=split/\t/,$gene;
	next if($strand ne $tmp[0]);
	if ($start>=$tmp[1] && $end <=$tmp[2]) {
		$ret .="\tintron";
		return $ret;
	}
		}
		$ret .="\tintergenic";
		return $ret;
	}
	return $ret;
}

sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -p -g -list -fa -species 
options:
-i input file,predictions file
-p input file,precursor sequence file
-g GeneListChrom input file
-list output file miRNA list file
-fa output file ,miRNA sequence fasta file.
-species  eg: hsa, ssc etc.
-h help
USAGE
exit(1);
}

