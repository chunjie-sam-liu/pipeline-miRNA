#!/usr/bin/perl -w
#Filename:
#Author: Jing Gong 
#Email:
#Date: 2011-12-27;
#Modified:
#Description:
my $version=1.00;

#use strict;
use Getopt::Long;
use File::Basename;


my %opts;
GetOptions(\%opts,"script=s","i=s","m=s","u=s","o=s","h");
if (!(defined $opts{script} and defined $opts{i} and defined $opts{m} and defined  $opts{u}  and defined $opts{o}  )|| defined $opts{h}) {
        &usage;
}

my $perlpath=$opts{'script'};
unless ($perlpath=~/\/$/){$perlpath.="/";}
my $mir_re=$opts{'i'};
my $mat_txt=$opts{'m'};
my $utr3=$opts{'u'};
my $outdir=$opts{'o'};
unless ($outdir=~/\/$/){$outdir.="/";}
my $sh=$outdir."target_predict.sh";

#task start

open OUTSH,">$sh";
print OUTSH "#!/bin/sh\n";

my $path = `echo \$PATH`;
print OUTSH "export PATH=$path\n\n";

open(IN,"<$mir_re") or die "cannot open $mir_re for target\n";
my ($mat,%hashmat);
while(<IN>){
	chomp;
	$_=~/(\S+),(\S+)/;
	$mat=$2;
	$hashmat{$mat}="";
}
close IN;

open(IN,"<$mat_txt") or die "cannot open $mir_re for target\n";
open(OUT,">$outdir"."tmp_seed.fa") ;
open(OUT2,">$outdir"."tmp_mat.fa");
while(<IN>){
	chomp;
my	@temp=split/\t/;
	if(exists $hashmat{$temp[0]}){
my 		$seed =substr($temp[1],1,7);
		print OUT "$temp[0]\t$seed\t9606\n"  ;         #here, 9606 is the human taxon; if other species, please change ;
		print OUT2 ">$temp[0]\n$temp[1]\n";
	}
}
close OUT;
close OUT2;
close IN;

open(IN,"<$utr3") or die "cannot open $utr3 for target\n";
open(OUT,">$outdir"."tmp_utr_for_targetscan") ;
while(<IN>){
	if(/>(\S+)/){
		
 print OUT "$1\t9606\t";
	$seq="";
	}
	while(<IN>){
	chomp;
	if(/>(\S+)/){
		
		print OUT "$seq\n$1\t9606\t";
	$seq="";
	}
	else{
	
		$seq.=$_;
		
	}
	}
}
print OUT "$seq\n";
close IN;
close OUT;


my $miranda_targets=$outdir."miranda_targets";

print OUTSH "miranda $outdir"."tmp_mat.fa $utr3 -sc 150 -en 15 -out $miranda_targets\n";
wait;

my $targetscan_targets=$outdir."TargetScan_targets";
print OUTSH "perl $perlpath"."targetscan_50.pl $outdir"."tmp_seed.fa $outdir"."tmp_utr_for_targetscan  $targetscan_targets\n";

wait;
my $target_predict=$outdir."target_predict.xls";
print OUTSH "perl $perlpath"."target_check_edit2.pl -mi $miranda_targets -ts $targetscan_targets  -o $target_predict\n\n";
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

