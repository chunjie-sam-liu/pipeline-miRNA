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

my ($infile)=@ARGV;
my (%opts, $outdir);
open (OPT,"<$infile") or die "cannot open config.properties\n";


while(my $line=<OPT>){
        if($line=~/RESULT_DIRECTORY\s*=\s*(\S+)/) { $opts{'out'}=$1; } 
		if($line=~/RAW_Data_PATH\s*=\s*(\S+)/)		{ $opts{'in'}=$1;}
	}
	close OPT;
	$data_path=dirname ( $opts{'in'});
	$opts{'p'}=$data_path."\/hsa_precursor.txt";                
	$opts{'m'}=$data_path."\/hsa_miRNA.txt"; 
	$opts{'fa'}=$opts{'out'}."\/ncrnamapping/mirna_candidate.fa";		
	$expoutfile=$opts{'out'}."\/miREXpress\/";

if (!(defined $opts{p} and defined $opts{m} and defined $opts{'fa'} ) ) { #necessary arguments
&usage;
}


my $pre=$opts{'p'};
my $mat=$opts{'m'};
my $fatxt=$opts{'out'}."\/ncrnamapping\/mirna_candidate.txt";
open OUT,">$fatxt";
open IN,"<$opts{'fa'}" or die "cannot open $opts{'fa'}"; #input file  
while (my $aline=<IN>) {
	chomp $aline;
	if($aline=~/>.*_x(\d+)/){
		print OUT "$1\t";
	}
	else{
		print OUT "$aline\n";
	}
}
if( ! -e $expoutfile){
mkdir $expoutfile;
}
my $sh=$opts{'out'}."\/miREXpress\/miRExpress.sh";
open OUTSH,">$sh" or die "cannot open $sh\n";
print OUTSH "/home/gongj/bin/miRExpress/src/alignmentSIMD -r $opts{'p'} -i $fatxt  -o $expoutfile \n\n";

print OUTSH "/home/gongj/bin/miRExpress/src/analysis -r  $opts{'p'} -m $opts{'m'} -d $expoutfile -o Example_alignment_result -t Example_expression_profile\n";

$var=$expoutfile."variation";
$nohit=$expoutfile."nohit";
print OUTSH "mkdir $var";
print OUTSH "/home/gongj/bin/miRExpress/src/alignmentSIMD -r  $opts{'p'} -i  $nohit -t 0.95 -o $var \n";
print OUTSH "/home/gongj/bin/miRExpress/src/analysis -r  $opts{'p'} -m $opts{'m'} -d $var -o Example_0.95_alignment_result -t Example_0.95_expression_profile\n";

`sh $sh`;
