#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:2010-11-09
#Description: É¾³ýmatched reads 
my $version=1.00;

#use strict;
use Getopt::Long;
use File::Basename;

my ($infile)=@ARGV;
my (%opts, $outdir);
open (OPT,"<$infile") or die "cannot open config.properties\n";


while(my $line=<OPT>){
	if($line=~/RESULT_DIRECTORY\s*=\s*(\S+)/)     { $opts{'out'}=$1;}
	if($line=~/GENOME_REFERENCE_FILE\s*=\s*(\S+)/){ $opts{'ref'}=$1;}
	if($line=~/MIRNA_REFERENCE_FILE\s*=\s*(\S+)/) { $opts{'micRNA'}=$1; }
	if($line=~/GENE_LIST_FILE\s*=\s*(\S+)/)       { $opts{'genelist'}=$1;}
	if($line=~/SPECIES_ALIAS_NAME\s*=\s*(\S+)/)   { $opts{'species'}=$1;}
	if($line=~/QUEUE_NAME\s*=\s*(\S+)/)           { $opts{'queue'}=$1;}
	if($line=~/PERL_PATH\s*=\s*(\S+)/)			  { $opts{'p'}=$1;}
}
close OPT;
		$opts{'o'}=$opts{'out'}."\/prediction\/";
		$opts{'nohit'}=$opts{'out'}."\/miREXpress\/nohit"; 
		$opts{'bp'}=$opts{'out'}."\/ncrnamapping\/mirna_candidate\.blastparsed"; 


if (!(defined $opts{'p'}) ) { #necessary arguments
&usage;
}

open(IN,"<$opts{'nohit'}") or die "cannot open nohit file\n";
$opts{'fa'}=$opts{'nohit'}."\.fa";
open(OUT,">$opts{'fa'}") or die "cannot open nohit.fa file\n";
while(<IN>){
chomp;

my($num,$seq)=split;
print OUT ">nohit",$.,"_x$num\n$seq\n";
}
close OUT;

$outdir=$opts{'o'};
unless ($outdir=~/\/$/) {$outdir .="/";}

if( ! -e $outdir){
`mkdir "$outdir"`;
}
my $perlpath=$opts{'p'};
#unless ($perlpath=~/\/$/){$perlpath .="/";}
my $queue = $opts{'queue'};

my $sh=$outdir."prediction.sh";

open OUTSH,">$sh" or die "cannot open $sh\n";
print OUTSH "#!/bin/sh\n";

my $path = `echo \$PATH`;
print OUTSH "export PATH=$path\n\n";

if(! -e $opts{'out'}."\/ncrnamapping\/ncrnamapfinish\.log"){
	my $sh2=$opts{'out'}."\/ncrnamapping\/ncrnamapping\.sh";
	print OUTSH "sh $sh2";
}
my $precursors=$outdir."precursors.fa";
print OUTSH "perl $perlpath"."excise_candidate.pl $opts{ref} $opts{bp} > $precursors\n\n";

my $structures=$outdir."structures";
print OUTSH "cat $precursors |RNAfold -noPS > $structures\n\n";

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
$0 config.properties 
USAGE
exit(1);
}

