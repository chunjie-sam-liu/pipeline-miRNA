#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:2010-11-09
#Description: É¾³ýmatched reads 
my $version=1.00;

use strict;
use Getopt::Long;
use File::Basename;


my %opts;
GetOptions(\%opts,"s=s","fa=s","bp=s","ref=s","micRNA=s","o=s","species=s","genelist=s","queue=s","h");
if (!(defined $opts{s} and defined $opts{fa} and defined $opts{bp} and defined $opts{o} and defined $opts{micRNA} and defined $opts{species} and defined $opts{genelist} and defined $opts{queue}) || defined $opts{h}) { #necessary arguments
&usage;
}

my $outdir=$opts{'o'};
unless ($outdir=~/\/$/) {$outdir .="/";}

my $perlpath=$opts{'s'};
unless ($perlpath=~/\/$/){$perlpath .="/";}
my $queue = $opts{'queue'};

my $sh=$outdir."prediction.sh";
my $qout=$outdir."stdout.txt";
my $qerr=$outdir."stderr.txt";

open OUTSH,">$sh";
print OUTSH "#!/bin/sh\n";
print OUTSH "#PBS -N mirna_prediction.sh\n";
print OUTSH "#PBS -o $qout\n";
print OUTSH "#PBS -e $qerr\n";
print OUTSH "#PBS -q $queue\n\n";

my $path = `echo \$PATH`;
print OUTSH "export PATH=$path\n\n";

my $precursors=$outdir."precursors.fa";
print OUTSH "perl $perlpath"."excise_candidate.pl $opts{ref} $opts{bp} > $precursors\n\n";
wait;

my $structures=$outdir."structures";
print OUTSH "cat $precursors |RNAfold --noPS > $structures\n\n";
wait;

my $tmp = $outdir."tmp";
my $signatures=$outdir."signatures";
print OUTSH "perl $perlpath"."auto_blast.pl $opts{fa} $precursors $tmp $perlpath -b > $signatures\n\n";
wait;

my $predictions=$outdir."predictions";
print OUTSH "perl $perlpath"."miRDeep.pl $signatures $structures $tmp -s $opts{micRNA} -y > $predictions\n\n";
wait;


print OUTSH "perl $perlpath"."miRDeep-conserve-novel.pl $opts{species} $predictions $opts{fa} $opts{bp} $opts{genelist} $outdir"."conserved.xls $outdir"."novel.xls $outdir"."remained.fasta $outdir"."loop.fasta $outdir"."mirna_prediction.stat $outdir"."miraln.result $outdir"."mirspe.result $outdir"."mirna.gff3\n\n";
wait;

my $finish = $outdir."predictfinish.log";
print OUTSH "touch $finish\n\n";


close OUTSH;

#system ("qsub $sh");
system ("sh $sh");
wait;

while(!(-e $finish)){
        sleep 60;
}



#task success finish
print "Finished running prediction jobs!";



sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -s -fa -bp -ref -o -micRNA -species
options:
-s script path
-fa #input file (reads align to no_mic_RNA in fasta format)
-bp #input file (reads align to mo_mic_RNA in BlastParsed format)
-ref #genome reference 
-micRNA  #fasta file with known miRNAs
-o #output directory
-species  #species name 
-queue #queue name
-h help
USAGE
exit(1);
}

