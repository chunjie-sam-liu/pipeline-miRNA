#!usr/bin/perl
#program:RPM_calculat.pl
#Author:zhanghm
#Date: Sat Jan 26 10:18:37 2013

#input reads file
#such as
#hsa-miR-576-3p  76
#hsa-miR-140-5p  197

use strict;
use Getopt::Long;
use File::Basename;

my($reads,$seq_depth,$rpm,$help);
GetOptions(
	'rpm=s' => \$rpm,
	'seq_depth=i' => \$seq_depth,
	'reads=s' => \$reads,
	'help|h' => \$help,
);

sub usage{
print <<USAGE;
usage:
#        perl $0 [options]
author
#        zhanghm   zhm.2009.happy@163.com
description:
		
options:
-h	--help	:print help info
-reads :(str): input reads number file
-seq_depth :(int): input the total reads number that map to genome
-rpm	:(str): output miRNA rpm file
e.g.:
#	perl $0 
USAGE
}

if(!(defined $reads && defined $rpm && defined $seq_depth) || defined $help)
{
	&usage();
	exit 0;
}

#read in miRNA reads file and build hash for them
open(IN,"<$reads") or die "can not open $reads\n";
open(OUT,">$rpm") or die "can not open $rpm\n";
print OUT "miRNA id\treads num\tRPM\n";
while(<IN>)
{
	chomp;
	my ($id,$num)=split("\t",$_);
	my $rpm = (1000000*$num)/$seq_depth;
	print OUT "$id\t$num\t";
	printf OUT "%.2f\n",$rpm;
}
close IN;
close OUT;
