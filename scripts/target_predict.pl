#!/usr/bin/perl -w
#Filename:
#Author:
#Email:
Date:
#Modified:
#Description:
my $version=1.00;

use strict;
use Getopt::Long;
use File::Basename;


my %opts;
GetOptions(\%opts,"script=s","i=s","p=s","g=s","species=s","e=s","u=s","s=s","o=s","que=s","h");
if (!(defined $opts{script} and defined $opts{i} and defined $opts{p} and defined $opts{g} and defined $opts{species} and defined $opts{e} and defined $opts{u} and defined $opts{s} and defined $opts{o} and defined $opts{que} )|| defined $opts{h}) {
        &usage;
}

my $perlpath=$opts{'script'};
unless ($perlpath=~/\/$/){$perlpath.="/";}
my $fa=$opts{'i'};
my $precuror=$opts{'p'};
my $genelist=$opts{'g'};
my $species=$opts{'species'};
my $energy=$opts{'e'};
my $utr3=$opts{'u'};
my $RNAhybrid_utr3=$opts{'s'};
my $outdir=$opts{'o'};
unless ($outdir=~/\/$/){$outdir.="/";}
my $sh=$outdir."target_predict.sh";
my $qout=$outdir."stdout.txt";
my $qerr=$outdir."stderr.txt";

#my $queue=$opts{'que'};

my $t_que=$opts{'que'};
my $q_que=$opts{'que'};



#task start

open OUTSH,">$sh";
print OUTSH "#!/bin/sh\n";
print OUTSH "#PBS -N target_prediction.sh\n";
print OUTSH "#PBS -o $qout\n";
print OUTSH "#PBS -e $qerr\n";
print OUTSH "#PBS -q $t_que\n\n";

my $path = `echo \$PATH`;
print OUTSH "export PATH=$path\n\n";

my $mirlist=$outdir."microRNA.list";
my $mirfa=$outdir."microRNA.fa";
print OUTSH "perl $perlpath"."miRNA_Express_and_sequence.pl -i $fa -p $precuror -g $genelist -list $mirlist  -fa $mirfa -species $species\n";

my $miranda_targets_out=$outdir."miranda_targets.out";
my $miranda_targets=$outdir."miranda_targets";
print OUTSH "miranda $mirfa $utr3 -sc 150 -en $energy -out $miranda_targets_out\n";
wait;
print OUTSH "perl $perlpath"."mi_result.pl -i $miranda_targets_out -o $miranda_targets\n";
wait;
my $RNAhybrid_targets_out=$outdir."RNAhybrid_targets.out";
my $RNAhybrid_targets=$outdir."RNAhybrid_targets";
print OUTSH "RNAhybrid -f 2,8 -c -s $RNAhybrid_utr3 -v 3 -u 3 -e $energy -q $mirfa -t $utr3 > $RNAhybrid_targets_out\n";
wait;
print OUTSH "perl $perlpath"."hy_result.pl -i $RNAhybrid_targets_out -o $RNAhybrid_targets\n";
wait;
my $target_predict=$outdir."target_predict";
print OUTSH "perl $perlpath"."target_check.pl -mi $miranda_targets -hy $RNAhybrid_targets -o $target_predict\n\n";
wait;

my $finish=$outdir."finish";
print OUTSH "touch $finish\n\n";

close OUTSH;
#system ("qsub $sh");
system ("sh $sh");
wait;


#task success finish


sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -script -i -p -g -species -e -u -s -o -que -h
options:
-script script path
-i	predictions file
-p      precurors file
-g      gene list chrom file
-species species alias name
-e	energy value
-u	utr3 file
-s      RNAhybrid_utr3
-o	out put directory
-que 	pbs queue name
-h	help
USAGE
exit(1);
}

