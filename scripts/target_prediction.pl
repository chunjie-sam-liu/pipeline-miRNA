#!/usr/bin/perl
#Author: liaoyifang 2012-6-18
#lyf890317@126.com
#script for finding mirna target gene

use Getopt::Long;
use File::Basename;

my %opts;
GetOptions(\%opts,"i=s","o=s","t=s","h");
#if(!(defined $opts{i} and defined $opts{o} and defined $opts{t}|| defined $opts{h})){
#	&usage;
#}

my $filein = $opts{'i'};
my $targetfile = $opts{'t'};
my $fileout = $opts{'o'};

open IN1, "<$filein"; #input file
open IN2, "<$targetfile"; #target file
open OUT, ">$fileout"; #output file

my @target = <IN2>;

while($line = <IN1>){
	chomp $line;
	if($line =~ /.*,(.*)\t(\d+)/){
		$mir1 = $1;
		$reads = $2;
		if($reads > 3){
			foreach $target(@target){
				chomp $target;
				if($target =~ /(.*)\t(.*)/){
					$mir2 = $1;
					$gene = $2;
					if($mir1 eq $mir2){
						print OUT $target."\n";
					}
				}
			}
		}
	}
}

close IN1;
close IN2;
close OUT;

#sub usage{
#	print <<"USAGE";
#	Usage:
#	perl target_prediction.pl -i infile -o outfile -t target_file
#	options:
#	-i input file
#	-o output file
#	-t target file
#	-h help
#	USAGE
#	exit(1);
#}
