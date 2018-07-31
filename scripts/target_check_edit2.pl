#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:
#Description: É¾³ýmatched reads 
my $version=1.00;

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"mi=s","ts=s","o=s","h");
if (!(defined $opts{mi} and defined $opts{ts} and defined $opts{o} ) || defined $opts{h}) { #necessary arguments
&usage;
}

my $filein=$opts{'mi'};
my $filein2=$opts{'ts'};
my $fileout=$opts{'o'};

#my %uniq;

my %hash;
my %hash2loc;
open IN,"<$filein"; #input file  
while (my $aline=<IN>) {
	
	if($aline!~/>/){ next;}
	if($aline=~m/>>/){ next;}
	chomp $aline;
	my @tmp=split/\t/,$aline;
	$tmp[0]=~s/>//;
	$hash{$tmp[0]}->{$tmp[1]}="$tmp[2]\t$tmp[3]\t$tmp[5]";
	 $tmp[5]=~/\d+\s+(\d+)/;
	$hash2loc{$tmp[0]}->{$tmp[1]}=$1;
	#print "$hash2loc{$tmp[0]}->{$tmp[1]}\n";
}
close IN;

my %target;
open OUT,">$fileout";
open IN ,"<$filein2";
while (my $aline=<IN>) {
	chomp $aline;
		my @tmp=split/\t/,$aline;	
		if (defined $hash{$tmp[1]}->{$tmp[0]} ) {
		my $loc2=abs(($hash2loc{$tmp[1]}->{$tmp[0]})-$tmp[6]);
		if($loc2<3){
		print OUT "$tmp[1]\t$tmp[0]\ttargetscan\t$tmp[5]\t$tmp[6]\tmiranda\t$hash{$tmp[1]}->{$tmp[0]}\r\n";
		}
		}	
	}
close IN;
close OUT;

sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -mi -hy -o
options:
-mi input file,miRanda target predict result
-hy input file,RNAhybrid target predict result.
-o output file
-h help
USAGE
exit(1);
}

