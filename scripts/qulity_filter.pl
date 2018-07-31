#!/usr/bin/perl
#12.07.09
#Author: zhuerle@163.com; edit by gongjing 2012-1-6;
#script for moving low quality and 3' 5' adapter and polyA

use strict;
use Getopt::Std;
use vars qw($opt_i $opt_o $opt_h);
getopts('i:o:h');
my $fq_file        = $opt_i;
my $outfile		   = $opt_o;
my $help			=$opt_h;

my $usage = << "USAGE";
Description: Perl script used to filter low quality short reads, remove polyA and trim 3' 5' adapter
Author: zhuerle\@163.com
Usage: perl Adapter_trim.pl [options] >outputfile
Options:
  -i <file>  Short reads file in fastq format 
  -o 		 outfile
Examples: perl Adapter_trim.pl -i sample.fq  -o outputfile
          perl Adapter_trim.pl -i sample.fq  -o outputfile
USAGE

if ($help) 
{
        print $usage;
		exit;
}

## fliter low quality
open( IN,"$fq_file") or die $!;
open (OUTFA,">$outfile" )or die $!;

while (my $line=<IN>)
{
	
	my $read = <IN>;
	my $third_line=<IN>;
	my $quality_line = <IN>;
	chomp $quality_line;
	 my $quality=0; 
	
	if(&check_qlt($quality_line,33,20)){
	print OUTFA $line.$read.$third_line.$quality_line."\n";
	}
}
close IN;
close OUTFA;


#########################################################################################################
#########################################################################################################
sub check_qlt
{
	my $quality_line = shift;
	my $asc = shift;
	my $tv = shift;
	my $num = 0;
	my $count = 0;
	my @ql = split (//,$quality_line);
	my $wid = $#ql+1;
	foreach my $i (0..$#ql)
	{
		$num = ord($ql[$i])-$asc;
		if($num-$tv <0){$count++;}
	}

		if( $count > $wid/2 ){return 0;}
		else {return 1;}
	
}


