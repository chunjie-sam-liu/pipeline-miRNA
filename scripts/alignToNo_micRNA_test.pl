#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:
#Description: É¾³ýmatched reads 
my $version=1.00;

#use strict;
use Getopt::Long;
use File::Basename;

my %opts;
GetOptions(\%opts,"s=s","bp=s","bf=s","nomiRNA=s","lf=s","o=s","queue=s","h");
if (!(defined $opts{s} and defined $opts{bp} and defined $opts{bf} and defined $opts{lf} and defined $opts{o} and defined $opts{s}  and defined $opts{nomiRNA} and defined $opts{queue}) || defined $opts{h})
{ #necessary arguments
&usage;
}


my $outdir=$opts{o};
unless ($outdir=~/\/$/) {$outdir.="/";}

my $log=$outdir."formatdb.log";
my $bp=$opts{'bp'};
my $fa=$opts{'bf'};
my $rfamlengthfile=$opts{'lf'};
my $sh=$outdir."ncrnamapping.sh";
my $qout=$outdir."stdout.txt";
my $qerr=$outdir."stderr.txt";
my $perlpath=$opts{'s'};
unless ($perlpath=~/\/$/){$perlpath .="/";} 
my $pos=5;
if(defined $opts{pos}){
	$pos=$opts{pos};
}
my $queue=$opts{'queue'};
open OUTSH,">$sh";
print OUTSH "#!/bin/sh\n";
print OUTSH "#PBS -N mirna_ncrnamapping.sh\n";
print OUTSH "#PBS -o $qout\n";
print OUTSH "#PBS -e $qerr\n";
print OUTSH "#PBS -q $queue\n\n";

my $path = `echo \$PATH`;
print OUTSH "export PATH=$path\n\n";
print OUTSH "formatdb -i $opts{nomiRNA} -p F -l $log\n";
wait;

my $blastout=$outdir."alignedToncrna.blastparsed";
print OUTSH "blastall -a 10 -p blastn -d $opts{nomiRNA} -i $fa -o $blastout -W 12 -D 2 -F F -m 8\n\n";
wait;

my $blast_ncrna_fa = $outdir."alignedToncrna.fa";
my $blast_ncrna_id = $outdir."alignedToncrna.id";
my $remain = $outdir."mirna_candidate";
my $count = $outdir."ncrna.count";

print OUTSH "perl $perlpath"."blastout_edit.pl -i $blastout -mis 1 -fa $fa -b $bp -id $blast_ncrna_fa -o $blast_ncrna_id -remain $remain >$count\n\n";
wait;

my $noMiRNA_count = $outdir."align_to_noMiRNA.count";
my $noMiRNA_txt= $outdir."no_miRNA.txt";
print OUTSH "perl $perlpath"."count_noMiRna_express.pl -i $blast_ncrna_id -l $rfamlengthfile -o $noMiRNA_count >$noMiRNA_txt\n\n";

my $mapsta = $outdir."counter_map.stat";
my $remain_bp = $remain.".blastparsed";

print OUTSH "perl $perlpath"."count_ncrna_map_reads.pl $bp $blast_ncrna_id $remain_bp $mapsta\n\n";
my $finish = $outdir."ncrnamapfinish.log";
print OUTSH "touch $finish\n";


close OUTSH;

#system ("qsub $sh");
system ("sh $sh");
wait;


while(!( -e $finish)){
	sleep 60;
}

sleep 60;

#task success finish
print "Finished running all jobs!";


sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -sam -genome -nomiRNA -o 
options:
-s script path
-sam  # .sam file 
-genome  genome reference
-nomiRNA no_mic_RNA sequence file
-o output directory
-queue queue
-h help
USAGE
exit(1);
}

