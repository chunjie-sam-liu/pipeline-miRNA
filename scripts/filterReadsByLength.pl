#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2010-01
#Modified:
#Description:  
my $version=1.00;

use strict;
use Getopt::Long;
use File::Basename;

my %opts;
GetOptions(\%opts,"i=s","tag:s","min=i","max=i","o=s","s=s","h");
if (!(defined $opts{i} and defined $opts{o} and defined $opts{min} and defined $opts{max}) || defined $opts{h}) { #necessary arguments
&usage;
}
open OUT,">$opts{o}";
open OUTSTA,">$opts{s}";
open IN,"<$opts{i}";
my $tag="fa";
$tag=$opts{'tag'} if(defined $opts{'tag'});

my (%uniquehash,%totalhash,@id,$num);
while (my $aline=<IN>) {
	chomp $aline;
	@id = split(/x/, $aline);
	$num = $id[1];
	my $seq=<IN>;
	chomp $seq;
	$uniquehash{length($seq)}++;
	$totalhash{length($seq)}+=$num;
	if (length ($seq)>=$opts{'min'}  && length ($seq) <=$opts{'max'}) {
		print OUT "$aline\n$seq\n";
	}
}
close IN;
close OUT;
#close OUTSTA;

print OUTSTA "#length\tunique\ttotal\n";
foreach  (sort{$a<=>$b} keys %uniquehash) {
	next if($_<$opts{'min'} || $_ > $opts{'max'});
	if($tag eq "fastq"){
		print OUTSTA $_,"\t",$uniquehash{$_}/2,"\t",$totalhash{$_}/2,"\n";
		next;
	}
	print OUTSTA $_,"\t",$uniquehash{$_},"\t",$totalhash{$_},"\n";
}
close OUTSTA;
sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -o  -min -max
options:

-i input file
-o output file
-s ouput statics file
-tag # the reads format  : fa or fastq,default fa.
-min reads min length.
-max reads max length.

-h help
USAGE
exit(1);
}

