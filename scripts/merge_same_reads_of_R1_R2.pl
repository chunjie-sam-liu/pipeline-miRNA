#!usr/bin/perl
#program:merge_same_reads_of_R1_R2.pl
#Author:zhanghm
#Date: Thu Jul 18 09:46:41 2013

use strict;
use Getopt::Long;
use File::Basename;

my($read1,$read2,$comb,$help);
GetOptions(
	'help|h' => \$help,
	'r1:s' => \$read1,
	'r2:s' => \$read2,
	'o:s' => \$comb,
);

sub usage{
print <<USAGE;
usage:
#        perl $0 [options]
author
#        zhanghm   zhm.2009.happy@163.com
description:
		this perl will merge the same reads from R1 and R2 clean fasta file
options:
-h	--help	:print help info
-r1	:input read 1 clean fasta file
-r2	:input read 2 clean fasta file
-o	:output merged result
e.g.:
#	perl $0 
USAGE
}

if(!(defined $read1 && defined $read2 && defined $comb ) || defined $help)
{
	&usage();
	exit 0;
}

#input read 1 file
open(IN,"$read1") or die $!;
my %merged_reads;
my $num1;
while(<IN>)
{
	chomp;
	if(/>.*x([0-9]+)/) 
	{
		$num1=$1;
	}
	else
	{
		$merged_reads{$_}+=$num1;
	}
}
close IN;

#input read 2 file
open(IN2,"$read2") or die $1;
my $num2;
while(<IN2>)
{
	chomp;
	if(/>.*x([0-9]+)/) 
	{
		$num2=$1;
	}
	else
	{
		$merged_reads{$_}+=$num2;	
	}
}
close IN2;

#print out the result
open(OUT,">$comb") or die $!;
my $num=0;
foreach my $key (sort {$merged_reads{$b} <=> $merged_reads{$a}} keys %merged_reads)
{
	print OUT ">sample_".$num."_x".$merged_reads{$key},"\n";
	print OUT $key,"\n";
	$num++;
}
close OUT;
