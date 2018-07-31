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
	if($aline=~/>>/){print OUT $aline;}
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

