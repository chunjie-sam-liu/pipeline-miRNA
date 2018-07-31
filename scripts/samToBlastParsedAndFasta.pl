#!/usr/bin/perl -w
#Filename:
#Author:		Yu Caixia
#Email:			yucx@big.ac.cn
#Date:		  	2010-5-5
#Modified:		
#Description:		Solexa Fragment reads, sam file (bwa result) to blastparsed_excision file and to fasta file
# 比上负链的，mapping位置转成Reference反向互补链上的位置
# 比上负链的，反向互补reads转成原Reads
# Fragment, unmapped, 一行有10列, unique mapped, 一行有18列, random mapped, 一行有19列
# paired,   unmapped, 一行有10列, unique mapped, 一行有20列, random mapped, 一行有21列
my $version=1.00;

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"i=s","j=s","o=s","p=s","h");
if (!(defined $opts{i} and defined $opts{j} and defined $opts{o} and defined $opts{p}) || defined $opts{h}) {
	&usage;
}

my $filein=$opts{'i'};
my $fileintwo=$opts{'j'};
my $fileout=$opts{'o'};
my $fileouttwo=$opts{'p'};
#my $alignment=$opts{'c'};

open IN,"<$filein";		#input file. fragment sam file (bwa result)
open INTWO,"<$fileintwo";	#input file. reference file
open OUT,">$fileout";		#output file. blastparsed_excision file
open OUTTWO,">$fileouttwo";	#output file. fasta file

my %ref_length_hash;
my ($id,$sequence)=();

while (my $aline=<INTWO>){
	chomp ($aline);
	if ($aline=~/^>(\S+)(.*)/){
		$id=$1;
		$sequence="";
		while (my $temp=<INTWO>){
			chomp ($temp);
			if ($temp=~/^>(\S+)(.*)/){
				$ref_length_hash{$id}=length($sequence);
				$id=$1;
				$sequence="";
				next;
			}
			$sequence.=$temp;
		}
	}
}
$ref_length_hash{$id}=length($sequence);

while (my $aline=<IN>){
	my ($query_length,$query_alignment_range,$ref_length,$ref_alignment_start,$ref_alignment_end,$ref_alignment_range,$strand,$reads);
	chomp ($aline);
	next if ($aline=~/^@/);
	my @info=split/\t/,$aline;
	if ($info[1]==4){
		next;
	}elsif($info[1]== 0 || $info[1]==16){
		#my @temp=split/;/,$info[$#info]; #print "$#temp\n";
		#if ($#temp <= ($alignment-1)){
		$reads=$info[9];
		if ($info[5]=~/(\d+)M/){
			$query_length=$1;
		}
		$query_alignment_range="1..".$query_length;
		$ref_length=$ref_length_hash{$info[2]};
		my $Ox = sprintf("%1x",$info[1]);
		if ($Ox eq "0"){
			$strand="+";
			$ref_alignment_start=$info[3];
			$ref_alignment_end=$info[3]+$query_length-1;
		}elsif ($Ox eq "10"){
			$strand="-";
			##### 比上负链的，mapping位置转成Reference反向互补链上的位置
			$ref_alignment_start=$ref_length-($info[3]+$query_length-1)+1;
			$ref_alignment_end=$ref_length-$info[3]+1;
			##### 比上负链的，反向互补reads转成原Reads
			$reads=~tr/ATGCatgc/TACGtacg/;
			$reads=reverse($reads);
		}
		$ref_alignment_range=$ref_alignment_start."..".$ref_alignment_end;
		print OUT "$info[0]\t$query_length\t$query_alignment_range\t$info[2]\t$ref_length\t$ref_alignment_range\t0.002\t1.00\t40.1\t$strand\n";
		print OUTTWO ">$info[0]\n$reads\n";

		next if($info[-1]!~/XA:Z:/);
		$info[-1]=~s/XA:Z://;
		my @temp=split/;/,$info[$#info];
		foreach  (@temp) {
			my @aa=split/,/,$_;
			if ($aa[2]=~/(\d+)M/){
				$query_length=$1;
			}
			$query_alignment_range="1..".$query_length;
			$ref_length=$ref_length_hash{$aa[0]};
			$aa[1]=~/([+-])(\d+)/;
			if ($1 eq "+") {
				$strand="+";
				$ref_alignment_start=$2;
				$ref_alignment_end=$2+$query_length-1;
			}
			else{
				$strand="-";
				##### 比上负链的，mapping位置转成Reference反向互补链上的位置
				$ref_alignment_start=$ref_length-($2+$query_length-1)+1;
				$ref_alignment_end=$ref_length-$2+1;
				##### 比上负链的，反向互补reads转成原Reads
				$reads=~tr/ATGCatgc/TACGtacg/;
				$reads=reverse($reads);
			}
			$ref_alignment_range=$ref_alignment_start."..".$ref_alignment_end;
			print OUT "$info[0]\t$query_length\t$query_alignment_range\t$aa[0]\t$ref_length\t$ref_alignment_range\t0.002\t1.00\t40.1\t$strand\n";
		}
	}	
}

sub usage{
	print <<"USAGE";
Version $version
Usage:
	$0 -i <input file> -j <input file> -o <output file> -p <output file> -c <interger>
options:
	-i input file. fragment sam file (bwa result)
	-j input file. reference file
	-o output file. blastparsed file
	-p output file. fasta file
	-c interger. Trim any sequences that are mapping more times than designed
	-h help
USAGE
	exit(1);
}

#perl samToBlastParseAndFasta.pl -i fragment.sam -j reference.fa -o re_fastare_blastparsed -p re_fasta -c 5
