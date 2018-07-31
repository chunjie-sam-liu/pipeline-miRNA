#!/usr/bin/perl 
#Author: Gong Jing 2011-12-5; 
#Description: used to compare common and specific reads of 2 dataset;
=infile data format: filter/sequence_cluster.txt
TGAGCATGTAGACAAAGGTAACACTGAAG   289574
GAGCATGTAGACAAAGGTAACACTGAAG    127103

usage: $0 $resultdirs  -------
such as:  perl perl/compare_miRExpress_express.pl /home/gongj/MIR_Pipeline/s2_result/ /home/gongj/MIR_Pipeline/s3_result/ /home/gongj/MIR_Pipeline/SR90_result/ /home/gongj/MIR_Pipeline/SR91_result/ 
please save the result file in the same dirctory. first parameter is the dirctory, following parameters are result dirctory;

=cut
use Statistics::R;
use File::Basename;

my ($dirs)=@ARGV; 
@files=split(/:/,$dirs);
my $dir=dirname ($files[0]);

$dir=~s/\/$//;
my $outdir=$dir."\/compare\/";                          # save the result in the /compare_result/ dirctory;
if( ! -e $outdir){
mkdir  $outdir ; 
}
open(OUT,">$outdir"."clean_read_total_compare.xls");
open(OUT2,">$outdir"."clean_read_uniq_compare.xls");
print OUT "sample";
print OUT2 "sample";
foreach $file (@files){
	open(IN,"$file"."/filter/sequence_cluster.txt") or die $!;
	$file=~s/\/.*\///g;
	print OUT "\t$file";
	print OUT2 "\t$file";
	while(<IN>){
		chomp;
		($reads,$num)=split;
#	print "$file\t$reads\t$num\n";
		$hash{$file}{$reads}=$num;
		$total{$file}+=$num;
		$uniq{$file}++;
		#	$hash2{$reads}{$file}="";
	}
	close IN;
print "$uniq{$file}\t$total{$file}\n";
}
print OUT "\n";
	print OUT2 "\n";


for($i=0;$i<=$#files;$i++){
	$files[$i]=~s/\/.*\///g;
	print OUT "$files[$i]";
	print OUT2 "$files[$i]";
		for($j=0;$j<=$#files;$j++){
			$common=0;
			$common_n=0;
			if($i==$j){
				print OUT "\t$total{$files[$i]}";
				print OUT2 "\t$uniq{$files[$i]}";
			}
			elsif($i<$j){
				print OUT "\t-";
				print OUT2 "\t-";
			}
			else{
				foreach my $k (keys %{$hash{$files[$i]}}){
				
				if(exists $hash{$files[$j]}{$k}){
					
							$common++;
							if($hash{$files[$i]}{$k}>=$hash{$files[$j]}{$k}){
								$common_n+=$hash{$files[$j]}{$k};
							}
							else{
								$common_n+=$hash{$files[$i]}{$k};
							}
						}
					}
				
			print OUT "\t$common_n";
			print OUT2 "\t$common";
			}
		}
		 print OUT "\n";
		  print OUT2 "\n";
	  }
