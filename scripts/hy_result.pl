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
GetOptions(\%opts,"i=s","o=s","h");
if (!(defined $opts{i} and defined $opts{o} ) || defined $opts{h}) { #necessary arguments
&usage;
}

my $filein=$opts{'i'};
my $fileout=$opts{'o'};

open IN,"<$filein"; #input file  
open OUT,">$fileout"; #output file  
while (my $aline=<IN>) {
	chomp $aline;
	my @temp=split/:/,$aline;
	next if(@temp<11);
	print OUT ">$temp[2]\t$temp[0]

 target length: $temp[1]
 miRNA length: $temp[3]
 mfe: $temp[4] kkcal/mol
 p-value: $temp[5]
 position: $temp[6]

 target 5\' $temp[7] 3\'
           $temp[8]
           $temp[9]
 miRNA  3\' $temp[10] 5\'


";
}

close IN;
close OUT;
sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -o
options:
-i input file
-o output file
-h help
USAGE
exit(1);
}

