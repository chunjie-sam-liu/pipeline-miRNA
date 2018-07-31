#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:2010-11-10
#Description: É¾³ýmatched reads 
my $version=1.00;

use strict;
use Getopt::Long;
use File::Basename;


my %opts;
GetOptions(\%opts,"i=s","s=s","a=s","min:i","max:i","o=s","queue=s","h");
if (!(defined $opts{i} and defined $opts{s} and defined $opts{o} and defined $opts{queue}) || defined $opts{h}) { #necessary arguments
&usage;
}


my $outdir=dirname($opts{'o'});
my $queue = $opts{'queue'}; 
my $perlpath=$opts{'s'};
unless ($perlpath=~/\/$/){$perlpath .="/";}

my $min = 15;
my $max = 35;
if(defined $opts{'min'}){
	$min = $opts{'min'}; 
}
if(defined $opts{'max'}){
	$max = $opts{'max'};
}
my $sh=$outdir."/preprocess.sh";
my $qout=$outdir."/stdout.txt";
my $qerr=$outdir."/stderr.txt";
my $finishfile = $outdir."/filterfinish.log";
open OUTSH,">$sh";
print OUTSH "#!/bin/sh\n";
#print OUTSH "#PBS -N mirna_preprocess.sh\n";		#here, add "#" for we donot use PBS systerm;
#print OUTSH "#PBS -o $qout\n";				  		#here, add "#" for we donot use PBS systerm;
#print OUTSH "#PBS -e $qerr\n";						#here, add "#" for we donot use PBS systerm;
#print OUTSH "#PBS -q $queue\n\n";					#here, add "#" for we donot use PBS systerm;

my $path = `echo -n \$PATH`; 						# here, i add -n; 
print OUTSH "export PATH=$path:$perlpath\n\n";

#my $stat = $opts{'o'}.".stat";
#open OUTFILE, "> $stat";

#my $lines=`wc -l $opts{i}`;
#my @temp=split/\s+/,$lines;

#print OUTFILE "Total reads:",$temp[0]/4,"\n" ; 

#my $adaper=$opts{'o'}."_adaper";
my $fa=$opts{'o'}."\.fa";    						  #here, mv line #86 here;
if (defined $opts{a}){
	
	print "\nadaper seq is $opts{a}\n";
#	print OUTSH "$perlpath"."fastx_clipper -a $opts{a} -i $opts{i} -o $adaper\n";
	print OUTSH "perl $perlpath"."Adapter_trim.pl -x 'GUUCAGAGUUCUACAGUCCGACGAUC' -y 'TCGTATGCCGTCTTCTGCTTG' -i $opts{i}  -o $fa\n";   #here, add this line to remove 5_adapter and 3_adapter and remove length <18;
	wait;
}
#else {$adaper=$opts{i};}

#my $cluster=$opts{o}."_cluster";

#print OUTSH "perl $perlpath"."groupReads.pl -input $adaper -o $cluster -tag solexa -minNr 0 -qual min:0\n\n";
#wait;

#print OUTSH "perl $perlpath"."filterReadsByLength.pl -i $cluster".".fa -o $opts{o}".".fa -s $opts{o}_clean".".stat -tag fa -min $min -max $max \n\n";
print OUTSH "touch $finishfile\n";

wait;


close OUTSH;

wait;
#system ("qsub $sh");
system ("sh $sh");
wait;
sleep 60;

#my $fa=$opts{'o'}.".fa";

while(!(-e $fa)){
      sleep 30;
}


#if(defined $opts{a}){
#       $lines=`wc -l $adaper`;
#        @temp=split/\s+/,$lines;
#		print OUTFILE "Adaptor removed:",$temp[0]/4,"\n" ;

#}


#$cluster = $cluster.".fa";	
#$lines=`wc -l $cluster`;
#@temp=split/\s+/,$lines;
#print OUTFILE "Clustered:",$temp[0]/2,"\n" ;

#$lines=`wc -l $fa`;
#@temp=split/\s+/,$lines;
#print OUTFILE "Selected tags:",$temp[0]/2,"\n" ;
#close OUTFILE;

#task success finish
sleep(30);
system("touch $finishfile");

sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -o -a -max -min
options:
-i input file
-s path of script
-o output file prefix
-a adaper sequence
-max micRNA reads max length
-min micRNA reads min length
-h help
USAGE
exit(1);
}

