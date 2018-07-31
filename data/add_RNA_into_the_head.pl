#!usr/bin/env perl

die "input two files: fasta, out\n" if(@ARGV!=2);
($in,$out)=@ARGV;

open(IN,"$in") or die $!;
open(OUT,">$out") or die $!;
while(<IN>)
{
	chomp;
	if(/^>.*RNA/){print OUT "$_\n"; $seq=<IN>; print OUT $seq;}
	else{ $new_id=$_."_ncRNA"; print OUT $new_id,"\n"; $seq=<IN>;print OUT $seq;}
}
close IN;
close OUT;
