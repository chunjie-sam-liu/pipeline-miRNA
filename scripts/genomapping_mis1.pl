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
use File::Basename;


my %opts;
GetOptions(\%opts,"s=s","sam=s","sai=s","genome=s","queryseq=s","queue=s","h");
if (!(defined $opts{s} and defined $opts{sam} and defined $opts{sai} and defined $opts{genome} and defined $opts{queryseq} and defined $opts{queue}) || defined $opts{h}) { #necessary arguments
&usage;
}

my $perlpath=$opts{'s'};
unless ($perlpath=~/\/$/){$perlpath .="/";}
my $outdir=dirname($opts{'sam'});
my $saifile=$opts{'sai'};
my $samfile=$opts{'sam'};
my $refseq=$opts{'genome'};
my $queryseq=$opts{'queryseq'};
my $queue = $opts{'queue'};
my $sh=$outdir."/bwa_aln.sh";
my $qout=$outdir."/stdout.txt";
my $qerr=$outdir."/stderr.txt";
my $log=$outdir."/genomemapfinish.log";
my $statfile=$outdir."/map_to_genome.count";
my $blastfile=$outdir."/alignedTogenome.blastparsed";
my $fafile=$outdir."/alignedTogenome.fa";
my $refseqfai= $refseq.".fai";
my $bamfile = $outdir."/sequence.bam";
my $bamsort = $outdir."/sequence";


open OUTSH,">$sh";
print OUTSH "#!/bin/sh\n";
print OUTSH "#PBS -N mirna_bwa_aln.sh\n";
print OUTSH "#PBS -o $qout\n";
print OUTSH "#PBS -e $qerr\n";
print OUTSH "#PBS -q $queue\n\n";

#my $path = `echo \$PATH`; changed by zhanghm
#print OUTSH "export PATH=$path\n\n"; changed by zhanghm
print OUTSH "bwa aln -n 1 -t 10 $refseq  $queryseq   >$saifile\n\n";
print OUTSH "bwa samse -f $samfile  $refseq  $saifile   $queryseq\n\n";
print OUTSH "samtools import $refseqfai $samfile $bamfile\n\n"; #add by tang 20110120
print OUTSH "samtools sort $bamfile $bamsort\n\n";
print OUTSH "samtools index $bamfile\n\n";

print OUTSH "perl $perlpath"."count_genome_map_reads.pl -i $samfile -o $statfile\n\n";

print OUTSH "perl $perlpath"."samToBlastParsedAndFasta.pl -i $samfile -j $refseq -o $blastfile -p $fafile\n\n";
print OUTSH "touch $log\n";
close OUTSH;

#system ("qsub $sh");
system ("sh $sh");
wait;

while(!(-e $log) ){
	sleep 60;
}

#task success finish
print "Finished running all jobs!";


#GetOptions(\%opts,"i=s","a=s","min:i","max:i","o=s","pid=i","h")
sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -s -sam -sai -genome -queryseq -pid 
options:
-s perl script path
-sam  # .sam file 
-sai sai file
-genome  genome reference
-queryseq queryseq
-h help
USAGE
exit(1);
}
