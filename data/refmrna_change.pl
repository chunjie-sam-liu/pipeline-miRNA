#!usr/bin/env perl

die "input two files: refMrna.fa, out" if(@ARGV!=2);
($in,$out)=@ARGV;

#input refmRNA
#>NM_001037228 1
#tcatgccaaaaggcagcaaataagtgccttttcttcccttcagaatacat
#ggacaatccaaagctctattagtct
open(IN,"$in") or die $!;
open(OUT,">$out") or die $!;
while(<IN>)
{
	chomp;
	if(/>(.*?)\s/)
	{
		$id=$1;
		print OUT "\n>$id"."_mRNA\n";
	}
	else
	{
		$_=~tr/atcg/ATCG/;
		print OUT;
	}
}
close IN;
close OUT;
