#!/usr/bin/perl 
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
GetOptions(\%opts,"i=s","mis=i","fa=s","b=s","o=s","id=s","remain=s","h");
if (!(defined $opts{i} and defined $opts{mis} and defined $opts{b} and defined $opts{fa} and defined $opts{remain} and defined $opts{o} and defined $opts{id}) || defined $opts{h}) { #necessary arguments
&usage;
}
my (%hash,%reads,%reads2,%reads3,%hashmir,%out,@count);

open IN,"<$opts{fa}";
while (my $aline=<IN>) {
	chomp $aline;
	$aline=~s/^>//;
	my $seq=<IN>;
	chomp $seq;
	$reads{$aline}=$seq;
	$reads3{$aline}=$seq;
}
close IN;

open OUT,">$opts{id}";
open O,">$opts{o}";
open IN,"<$opts{i}";
my $tag="";
my $tag2="";
# sample_0_x289574        chromosome.1;processed_transcript;ENSG00000208317       100.00  29      0       0       1       29      35      63      8e-09   58.0

while (my $aline=<IN>) {
	chomp $aline;
	my @temp=split/\t/,$aline;
	if ($tag eq $temp[0]) {
		next;
	}
	$tag=$temp[0];
#	$tag2=$temp[11];
	my $mis=$temp[4]+length($reads{$temp[0]})-$temp[3];
	if ($mis<=$opts{'mis'} && $temp[8]<$temp[9]) {
		print OUT ">$temp[0]\n$reads{$temp[0]}\n";
		print O "$temp[0]\n>$temp[1]\n";
		$count[0]++;
		$temp[0]=~/_x(\d+)/;
		$count[1]+=$1;
	#	if($temp[1]=~/miR|mir|let/){
	#		$hashmir{$temp[0]}=$reads3{$temp[0]};
	#	}
		if($temp[1]=~/rRNA|tRNA|snRNA|snoRNA/ ){
		delete ($reads{$temp[0]});
		
	}
	}
}

close IN;
close OUT;
close O;
my %reads2=(%reads,%hashmir);
my $fa=$opts{'remain'}.".fa";
my $blast=$opts{'remain'}.".blastparsed";
open F,">$fa";
open B,">$blast";
open IN,"<$opts{b}";
while (my $aline=<IN>) {
	chomp $aline;
	my @temp=split/\t/,$aline;
	if (defined $reads2{$temp[0]}) {
		print B "$aline\n";
	}
}
close IN;
close B;
foreach my $key (keys %reads2) {
	print F ">",$key,"\n",$reads2{$key},"\n";
	$count[2]++;
	$key=~/_x(\d+)/;
	$count[3]+=$1;
}
close F;

print "#class\tno_miRNA\tremain\n";
print "#unique_map\t$count[0]\t";
print "$count[2]\n";
print "#total_map\t$count[1]\t";
print "$count[3]\n";
#print "#unique_map\t$count[0]\t$count[2]\n";
#print "#total_map\t$count[1]\t$count[3]\n";
sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -o -mis -fa -id -b -remain
options:
-i input file #blast out file (-m 8 )
-fa input file # fa file map to genome
-b input file # blastparsed file map to genome
-mis max number mismatch
-id output file # fa format file,reads align to nomiRNA.
-o output file #id_align_to_nomiRNA
-remain output file prefix # reads remained print out in .fa format and .blastparsed format.
-h help
USAGE
exit(1);
}

