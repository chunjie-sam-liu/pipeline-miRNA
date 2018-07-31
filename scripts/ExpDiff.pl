#!/usr/bin/perl

use strict;
use Cwd 'abs_path';
use Getopt::Long;

=head1 Description:
	ExpDiff.pl
	Find expression difference between genes

=head1 Date:
	2012-12-18

=head1 Usage:
	perl ExpDiff.pl [options]

	options:
	-s1: sample 1 miRNA expression
	-s2: sample 2 miRNA expression
	-log2: log2 cutoff of fold-change, default is 1
	-fdr: fdr cutoff, default is 0.01
	-tpm: tpm cutoff, default is 5
	-out: output file
	-help|?: help informaition

	e.g.: perl ExpDiff.pl -compare diff_config.list -outdir .

=cut

my %opts;
GetOptions(\%opts,"s1=s","s2=s","log2=f","fdr=f","tpm=f","outdir=s","perl_path=s","h");
if (!(defined $opts{s1} and defined $opts{s2} and defined $opts{outdir} and defined $opts{perl_path}) || defined $opts{h})
{
	&usage;
}

my ($log2,$fdr,$tpm);
$log2 = $opts{"log2"} ? $opts{"log2"} : 1;
$fdr = $opts{"fdr"} ? $opts{"fdr"} : 0.01;
$tpm = $opts{"tpm"} ? $opts{"tpm"} : 5 ;
#print $log2,"\t",$fdr,"\t",$tpm,"\n";

#sample expression profile file
#such as
#ID      Total tags      TMP
#ssc-miR-378     372     472.67
#ssc-miR-339     95      120.71
#ssc-miR-103     7626    9689.68

my %genes = (); ## len, readsnum, RPKM
my %sampleReadsNum = ();
$opts{"s1"}=~/(.*?)_/;
my $s1_name=$1;
$opts{"s2"}=~/(.*?)_/;
my $s2_name=$1;
#print $s1_name,"\t",$s2_name,"\n";
my $sample1=$opts{"s1"};
my $sample2=$opts{"s2"};

#read in sample1 miRNA expression file
open(IN,"$sample1") or die "can not open $sample1\n";
<IN>;
while(<IN>)
{
	chomp;
	my @data=split;
	$sampleReadsNum{"sample1"}+=$data[1];
	$genes{$data[0]}{"sample1"}{"reads"} = $data[1];
	$genes{$data[0]}{"sample1"}{"tpm"} = $data[2];
}
close IN;

#read in sample2 miRNA expression file
open(IN2,"$sample2") or die "can not open $sample2\n";
<IN2>;
while(<IN2>)
{
	chomp;
	my @data=split;
	$sampleReadsNum{"sample2"}+=$data[1];
	$genes{$data[0]}{"sample2"}{"reads"} = $data[1];
	#print $genes{$data[0]}{"sample2"}{"reads"},"\n";
	$genes{$data[0]}{"sample2"}{"tpm"} = $data[2];
}
close IN2;

# diff
my %diffs = ();
my ($total,$sort_i);  # total genes, order number of sort genes
if( $sampleReadsNum{"sample1"}<=0 || $sampleReadsNum{"sample2"}<=0 )
{
	die "total reads number of $s1_name or $s2_name is 0\n";
}

# p-value
my $tmp=$opts{"outdir"}."/S1_vs_S2";
open my $PV, ">$tmp" || die $!;
for my $gene (keys %genes) {
	next if ((!exists $genes{$gene}{"sample1"} || $genes{$gene}{"sample1"}{"reads"} == 0) && (!exists $genes{$gene}{"sample2"} || $genes{$gene}{"sample2"}{"reads"} == 0));
	if (!exists $genes{$gene}{"sample1"} or $genes{$gene}{"sample1"}{'reads'} == 0) {
		$genes{$gene}{"sample1"}{'reads'} = 0;
		$genes{$gene}{"sample1"}{'tpm'} = 0.001;
	}
	if (!exists $genes{$gene}{"sample2"} or $genes{$gene}{"sample2"}{'reads'} == 0) {
		$genes{$gene}{"sample2"}{'reads'} = 0;
		$genes{$gene}{"sample2"}{'tpm'} = 0.001;
	}
	print $PV "$gene\t$sampleReadsNum{'sample1'}\t$sampleReadsNum{'sample2'}\t$genes{$gene}{'sample1'}{'reads'}\t$genes{$gene}{'sample2'}{'reads'}\n";
}
close $PV;
system("$opts{perl_path}/statistic -i $tmp -o $opts{outdir}/GeneDiffExp.xls && rm $tmp");

open my $IN, "<$opts{outdir}/GeneDiffExp.xls" || die $!;
while (<$IN>) 
{
	chomp;
	my @tabs = split /\t/, $_;
	$diffs{$tabs[0]}{"pvalue"} = $tabs[5] * 2;
}
$total = $.;    #total genes number
close $IN;
system("rm $opts{outdir}/GeneDiffExp.xls");

# log2, fdr, up/down
$sort_i = 0;    #order number of gene
for my $gene (sort {$diffs{$a}{"pvalue"} <=> $diffs{$b}{"pvalue"}} keys %diffs) 
{
	$sort_i++;
	$diffs{$gene}{"fdr"} = $diffs{$gene}{"pvalue"} * $total / $sort_i;
	my $tpmA = $genes{$gene}{"sample1"}{"tpm"};
	my $tpmB = $genes{$gene}{"sample2"}{"tpm"};
	#print "$tpmA\t$tpmB\n";
	$diffs{$gene}{"log2"} = log($tpmB/$tpmA)/log(2);
	if ($diffs{$gene}{"log2"} > 0) 
	{
		$diffs{$gene}{"ud"} = "Up";
	}
	elsif ($diffs{$gene}{"log2"} < 0) 
	{
		$diffs{$gene}{"ud"} = "Down";
	}
	else 
	{
		$diffs{$gene}{"ud"} = "-";
	}
}

# print result
open my $OUT, "> $opts{outdir}/$s1_name"."_vs_".$s2_name."_GeneDiffExp.xls" || die $!;
open my $FT, "> $opts{outdir}/$s1_name"."_vs_".$s2_name."_GeneDiffExpFilter.xls" || die $!;

my $head = "geneID\t$s1_name-Expression\t$s2_name-Expression\t$s1_name-TPM\t$s2_name-TPM\tsubstraction\tlog2 Ratio($s2_name/$s1_name)\tUp-Down-Regulation($s2_name/$s1_name)\tP-value\tFDR\n";
print $OUT $head;
print $FT $head;

for my $gene (sort {sgn($diffs{$b}{"log2"}) <=> sgn($diffs{$a}{"log2"}) or abs($diffs{$b}{"log2"}) <=> abs($diffs{$a}{"log2"})} keys %diffs) 
{
	# sort, sort (1, 2, 0, -1, -2) to (2, 1, 0, -2, -1)
	my $diff=$genes{$gene}{sample2}{'tpm'}-$genes{$gene}{sample1}{'tpm'};
	#my $line = "$gene\t$genes{$gene}{sample1}{'reads'}\t$genes{$gene}{sample2}{'reads'}\t$genes{$gene}{sample1}{'tpm'}\t$genes{$gene}{sample2}{'tpm'}\t$diff\t$diffs{$gene}{'log2'}\t$diffs{$gene}{'ud'}\t$diffs{$gene}{'pvalue'}\t$diffs{$gene}{'fdr'}\n";
	printf $OUT "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t",$gene,$genes{$gene}{sample1}{'reads'},$genes{$gene}{sample2}{'reads'},$genes{$gene}{sample1}{'tpm'},$genes{$gene}{sample2}{'tpm'},$diff,$diffs{$gene}{'log2'},$diffs{$gene}{'ud'};
	print $OUT $diffs{$gene}{'pvalue'},"\t",$diffs{$gene}{'fdr'},"\n";
	#printf $OUT "%.2f\t%.2f\n",$diffs{$gene}{'pvalue'},$diffs{$gene}{'fdr'};
	# filter
	if (abs($diffs{$gene}{'log2'}) >= $log2 && $diffs{$gene}{'fdr'} <= $fdr && ($genes{$gene}{sample1}{'tpm'}>=$tpm || $genes{$gene}{sample2}{'tpm'}>=$tpm))
	{
		printf $FT "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t",$gene,$genes{$gene}{sample1}{'reads'},$genes{$gene}{sample2}{'reads'},$genes{$gene}{sample1}{'tpm'},$genes{$gene}{sample2}{'tpm'},$diff,$diffs{$gene}{'log2'},$diffs{$gene}{'ud'};
		print $FT $diffs{$gene}{'pvalue'},"\t",$diffs{$gene}{'fdr'},"\n";
	} 
		#	printf $FT "%.2f\t%.2f\n",$diffs{$gene}{'pvalue'},$diffs{$gene}{'fdr'};
}
close $OUT;
close $FT;

print STDOUT "Finish!\n";

############################################################################################
#substrine
sub sgn 
{
	my $n = shift;
	if ($n > 0)
	{
		$n = 1;
	}
	elsif ($n < 0)
	{
		$n = -1;
	} 
	else 
	{
		$n = 0;
	}
	return $n;
}

sub usage
{
print <<"USAGE";
Usage:
$0  -s1 -s2 -outdir
options:
-s1	input sample1
-s2	input sample2
-log2	input cutoff value for diffrential expression: default is 1
-fdr	input FDR value: default is 0.01
-tpm	input the minnum of TPM
-outdir	input comparation result
-perl_path	input perl path
-h	help
USAGE
exit(1);
}
