#!/usr/bin/perl -w
#Filename:
#Author: gongjing
#Email: doremi1985@sina.com
#Date: 2012-05-16
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
		if($line=~/PRE_MIRNA\s*=\s*(\S+)/)  {$opts{'p'}=$1;}
		if($line=~/MAT_MIRNA\s*=\s*(\S+)/)  {$opts{'m'}=$1;}
		if($line=~/3_ADAPTOR_FILE\s*=\s*(\S+)/) {$adaptor3=$1;}
		if($line=~/5_ADAPTOR_FILE\s*=\s*(\S+)/) {$adaptor5=$1;}

	}
	close OPT;
	$expoutfile=$opts{'out'}."\/miREXpress\/";

if (!(defined $opts{p} and defined $opts{m} and defined $opts{'in'} ) ) { #necessary arguments
&usage;
}


my $pre=$opts{'p'};
my $mat=$opts{'m'};

if( ! -e $opts{'out'}){
mkdir $opts{'out'};

mkdir  $expoutfile;
}
my $sh=$opts{'out'}."\/miRExpress.sh";
open OUTSH,">$sh" or die "cannot open $sh\n";
print OUTSH "perl /home/gongjing/MIR-pipeline/tools/perl/qulity_filter.pl -i  $opts{'in'} -o $opts{'in'}.qlt\n";
print OUTSH "/home/gongjing/tools/miRExpress/src/Raw_data_parse -r  $opts{'in'}.qlt -o  $opts{'in'}.qlt.merge\n";
print OUTSH "/home/gongjing/tools/miRExpress/src/Trim_adapter -i  $opts{'in'}.qlt.merge -t $adaptor3  -h $adaptor5  -o $opts{'in'}.qlt.merge.trim\n";
print OUTSH "/home/gongjing/tools/miRExpress/src/alignmentSIMD -r $opts{'p'} -i $opts{'in'}.qlt.merge.trim  -o $expoutfile \n\n";

print OUTSH "/home/gongjing/tools/miRExpress/src/analysis -r  $opts{'p'} -m $opts{'m'} -d $expoutfile -o Example_alignment_result -t Example_expression_profile\n";

$var=$expoutfile."variation";
$nohit=$expoutfile."nohit";
print OUTSH "mkdir $var";
print OUTSH "/home/gongjing/tools/miRExpress/src/alignmentSIMD -r  $opts{'p'} -i  $nohit -t 0.95 -o $var \n";
print OUTSH "/home/gongjing/tools/miRExpress/src/analysis -r  $opts{'p'} -m $opts{'m'} -d $var -o Example_0.95_alignment_result -t Example_0.95_expression_profile\n";

`sh $sh`;
