#!/usr/bin/perl 
#Author: Gong Jing 2011-12-5; 
#Description: used to compare miRNA expression;
=infile Example_expression_profile format
hsa-mir-1294,hsa-miR-1294       50

usage: $0  result_file1 result_file2 -------
such as:  perl perl/compare_miRExpress_express.pl /home/gongj/MIR_Pipeline/s2_result/ /home/gongj/MIR_Pipeline/s3_result/ /home/gongj/MIR_Pipeline/SR90_result/ /home/gongj/MIR_Pipeline/SR91_result/ 
please save the result file in the same dirctory. first parameter is the dirctory, following parameters are result dirctory;

=cut
use Statistics::R;
use File::Basename;


my @files=@ARGV; 
my $dir=dirname ($files[0]);

$dir=~s/\/$//;
my $outdir=$dir."\/compare\/";                          # save the result in the /compare_result/ dirctory;
if( ! -e $outdir){
mkdir  $outdir ; 
}

my $i=0;
my (%data,%mirs);
my ($infile,$infile2);
for($i=0;$i<=$#files;$i++){
	$files[$i]=~s/\/$//;
	$infile=$files[$i]."\/miREXpress\/Example_expression_profile";
	$infile2=$files[$i]."\/filter\/sequence.stat";
	open(IN,"<$infile") or die "cant open $infile\n"; 
	while (<IN>){
		chomp;
		$_=~/(\S+),(\S+)\s+(\d+)/;
	my	$pre=$1;
	my  $mir=$2;
	my  $exp=$3;
		if($exp>=3){
		$data[$i]{$mir}{'pre'}=$pre;
		$data[$i]{$mir}{'exp'}=$exp;
		$data2{$mir}{$i}=$exp;
		$mirs{$mir}="";	
		#	print "$data[$i]{$mir}{'chr'}\t $data[$i]{$mir}{'mat'}\t$data[$i]{$mir}{'star'}\n";
		}
	}
	close IN;
	open(IN2,"<$infile2") or die "cannot open $infile2\n";
	while(<IN2>){
		if($_=~/after\s+rm\s+polyA:\s*(\d+)/){$total[$i]=$1;}           # count the total number of each data set;
	}
	close IN2;	
#print "$total[$i]\n";
}
my	($total_c,$count,$total_c,$rpm);
my $outfile=$outdir."common_Known_miRNA_result\.xls"; 

open OUT,">$outfile" or die "cannot open $outfile\n"; 


print OUT "miRNA\t";
for($n=0;$n<=$#files;$n++){
	$file_name = basename $files[$n];
	print OUT "$file_name","_c\t$file_name","_rpm\t";
	$fh=OUTS.$n;
	$s=$n+1;
	open $fh,">$outdir"."sample_".$s."_specific_miRNAs.xls";
	print $fh "miRNAs\treads\ttpm\n";
	
	foreach my  $mir3 (sort keys %mirs){
	$k=keys %{$data2{$mir3}};
	if(exists $data2{$mir3}{$n} && $k==1 ){
	$total_c=$total[$n];
	$count=$data2{$mir3}{$n};
	$rpm =$count/$total_c*1000000;
	$rpm =sprintf ("%8.2f",$rpm);
	print $fh "$mir3\t$data2{$mir3}{$n}\t$rpm\n";
	}
}
}
print OUT "\r\n";
foreach my  $mir2 (sort keys %mirs){
	$k=keys %{$data2{$mir2}};
	if($k==($#files+1)){	
		print OUT "$mir2\t";
		for($j=0;$j<=$#files;$j++){
		$total_c=$total[$j];
		$count=$data2{$mir2}{$j};
		$rpm =$count/$total_c*1000000;
		 $rpm =sprintf ("%8.2f",$rpm);
	 print OUT "$count\t$rpm\t";
	 }

print OUT "\r\n";
}
}
close OUT;
wait;
## different expression miRNAs;
open(IN,"<$outfile") or die $!;
$outfile_diff=$outdir."different_expressed_miRNAs";
open(OUT,">$outfile_diff") or die $!;

while(<IN>){
	if(/miRNA/){ print OUT "$_"; next;}
	$tag=0;
	chomp;
my	@data=split/\t/;
	for($j2=0;$j2<=$#files;$j2++){
		$num=$data[$j2*2+2];
		for($j3=$j2+1;$j3<=$#files;$j3++){
			$num2=$data[$j3*2+2];
				if($num/$num2>=2 or $num/$num2<=0.5){
							$tag=1;
				}
		}
	}
	if($tag==1){
		print OUT  "$_\n";
	}
#	print "\n";
}
close IN;
wait;


####plot;
my $R = Statistics::R->new() ;
$R->startR ;
$com_pdf=$outdir."expression_heatmap.pdf";
$diff_pdf=$outdir."different_expressed_miRNA.pdf";
$s=$#files+1;
my $r=<<END;
n<-(1:$s)*2
colors<- colorRampPalette(c("white", "darkblue"))(100)
data.level1 <- read.table("$outfile", head=T, sep="\\t",row.names=1)
data2<-read.table("$outfile_diff", head=T, sep="\\t",row.names=1)
data.mc <- data.level1[,c(n)]
data.mc2 <- data2[,c(n)]
matrix.mc <- as.matrix(data.mc)
matrix.mc2 <- as.matrix(data.mc2)
pdf(file="$com_pdf")
heatmap(log2(matrix.mc), col=colors, scale="none", margins=c(5,9), cexCol=1.2)
graphics.off()\n
pdf(file="$diff_pdf")
heatmap(log2(matrix.mc2), col=colors, scale="none", margins=c(5,9), cexCol=1.2)
graphics.off()\n
END

$R->run($r);
$R->stop();

