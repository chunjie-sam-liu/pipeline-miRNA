#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2010-11-23
#Modified:
#Description:Compute no-miRNA express
my $version=1.00;

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"i=s","l=s","o=s","h");
if (!(defined $opts{i} and defined $opts{l} and defined $opts{o} ) || defined $opts{h}) { #necessary arguments
&usage;
}

my $filein=$opts{'i'};
my $filelength=$opts{'l'};
my $fileout=$opts{'o'};

#print $filein;
#print $filelength;
open IN,"<$filein"; #input file  
open INLEN,"<$filelength";#input file
open OUT,">$fileout"; #output file  
my %hash;my %hashlength;my @counter;my @total;
$counter[0]=$counter[1]=$counter[2]=$counter[3]=$counter[4]=$counter[5]=$counter[6]=$counter[7]=0;
$total[0]=$total[1]=$total[3]=$total[2]=$total[4]=$total[5]=$total[6]=$total[7]=0;
while (my $aline=<IN>) {
	chomp $aline;
	my $rna=<IN>;
	$aline=~/_x(\d+)$/;
	$hash{$rna}[0] ++;
	$hash{$rna}[1] +=$1;
}
close IN;

#fetch rna length
my @length;
while (my $lline=<INLEN>){
	chomp $lline;
	#print $lline;
	@length=split(/->/, $lline);
	$hashlength{$length[0]}=$length[1];
}
close INLEN;


print OUT "#Class\tRNA_Tag\tRNA_ID\tlength\tuniqe\ttotal\n";
foreach my $key (keys %hash) {
	$key=~/^>(\S+)/;
	#print $1."\n";
	my $temp=$1;
	print OUT $temp,"\t",$hashlength{$1},"\t",$hash{$key}[0],"\t",$hash{$key}[1],"\n";
	if ($temp =~ /mRNA/i){

		$counter[1]++; 
		$total[1]+=$hash{$key}[1];
	}elsif($temp =~ /tRNA/i || $temp =~ /transfer RNA/i){
		$counter[2]++;
		$total[2]+=$hash{$key}[1];
	}elsif($temp =~ /snoRNA/i ){	
		$counter[4]++;
		$total[4]+=$hash{$key}[1];
	}elsif($temp =~ /rRNA/i || $temp =~ /ribosomal RNA/i){
		$counter[3]++;
		$total[3]+=$hash{$key}[1];
	}elsif($temp =~ /snRNA/i){
		$counter[5]++;
		$total[5]+=$hash{$key}[1];
	}elsif($temp =~ /piRNA/i){
		$counter[6]++;
		$total[6]+=$hash{$key}[1];
	}
	elsif($temp =~/mir|let/i){
		$counter[7]++;
		$total[7]+=$hash{$key}[1];
	}
	elsif($temp =~ /lncRNA/i){
		$counter[8]++;
		$total[8]+=$hash{$key}[1];
	}
	else{
		$counter[9]++;
		$total[9]+=$hash{$key}[1];
	}
$total[0]+=$hash{$key}[1];

}

close OUT;

$counter[0] = $counter[1] + $counter[2] + $counter[3] + $counter[4] + $counter[5] + $counter[6] + $counter[7] + $counter[8]+$counter[9];
print "#the number of expressed RNA:\n";
print "mRNA:\t$counter[1]\t$total[1]\t (includes mRNA)\n";
print "tRNA:\t$counter[2]\t$total[2]\t (includes tRNA, transfer RNA)\n";
print "rRNA:\t$counter[3]\t$total[3]\t (includes ...rRNA..., ribosomal RNA, ...S(...))\n";
print "snoRNA:\t$counter[4]\t$total[4]\t (includes snoRNA )\n";
print "snRNA:\t$counter[5]\t$total[5]\t(includes snRNA)\n";
print "piRNA:\t$counter[6]\t$total[6]\t (includes piRNA)\n";
print "lncRNA:\t$counter[8]\t$total[8]\t (includes lncRNA)\n";
print "miRNA:\t$counter[7]\t$total[7]\t (includes miRNA)\n";
print "others:\t$counter[9]\t$total[9]\n";
print "total:\t$counter[0]\t$total[0]\n";

sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -o
options:
-i input file
-lengthfile input file
-o output file
-h help
USAGE
exit(1);
}

