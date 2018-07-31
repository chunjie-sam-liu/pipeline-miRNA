#!usr/bin/perl
#program:stat_miRNA_reads_by_alignment_result.pl
#Author:zhanghm
#Date: Mon Jul 22 15:30:07 2013

use strict;
use Getopt::Long;
use File::Basename;

my($align,$reads,$help);
GetOptions(
	'help|h' => \$help,
	'i:s' => \$align,
	'o:s' => \$reads,
);

sub usage{
print <<USAGE;
usage:
#        perl $0 [options]
author
#        zhanghm   zhm.2009.happy@163.com
description:
		this perl will calculate miRNA reads number by sample_alignment_result
options:
-h	--help	:print help info
-i	:input sample_alignment_result file
-o	:output miRNA reads count
e.g.:
#	perl $0 
USAGE
}

if(!(defined $align && defined $reads) || defined $help)
{
	&usage();
	exit 0;
}

#read in alignment file
#hsa-mir-548z
#AAGUAUUAAGUUGGUGCAAAAGUAAUUGAGAUUUUUGCUACUGAAAGUAAUGGCAAAAACCGCAAUUACUUUUGCACCAACCUAAUAGAUGCCAAUG
#                                                     CAAAAACCGCAAUUACUUUUGCA hsa-miR-548z
#                                                     CAAAAACCGCAATTACTTTTG 15
#                                                     CAAAAACCGCAATTACTTTTGC 11
#Others*******************************************************************************************
#                                                               AATTACTTTTGCACCAACCTAA 5
#//
open(IN,"<$align") or die $!;
$/="//";
my %hash;
while(<IN>)
{
	chomp;
	my $miRNA_id="space";
	my @data=split("\n",$_);
	for(my $i=0;$i<=$#data;$i++)
	{
		if($data[$i]=~/[ACGU]+\s+(.*miR.*|.*let.*)/)
		{
			$miRNA_id=$1;
			#print $miRNA_id,"\n";
		}
		elsif($data[$i]=~/[ACGT]+\s+\d+/)
		{
			$data[$i]=~s/^\s+//g;
			#print $data[$i],"\n";
			$hash{$miRNA_id}->{$data[$i]}=1;
		}
		elsif($data[$i]=~/Others/)
		{
			last;
		}
	}
}
close IN;
$/="\n";

#print out the result
open(OUT,">$reads") or die $!;
foreach my $key (keys %hash)
{
	my $num=0;
	foreach my $key2 (keys %{$hash{$key}})
	{
		$key2=~/[ACGT]+\s(\d+)/;
		$num+=$1;
		if($key eq "space")
		{
			print $key2,"\n";
		}
	}
	print OUT "$key\t$num\n";
}
close OUT;
