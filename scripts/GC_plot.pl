#!/usr/bin/perl
#12.07.09
#Author: zhuerle@163.com; edit by gongjing 2012-1-6;
#script for moving low quality and 3' 5' adapter and polyA

use Statistics::R;
use File::Basename;
use strict;

#num    seq
#9       ACACGATCACTCCCGCTGAGC
#1       CTGAGAAAATGGGATC

($infile,$outfile)=@ARGV;
open(IN,"<$infile") or die $!;
open(OUT,">$outfile") or die $!;
#%count="";
while(<IN>){
chomp;
 my($seq,$num)=(split/\t/)[1,2];
 @data=split(//,$seq);

 for($i=0;$i<=$#data;$i++){
 $count{$data[$i]}->{$i}=$count{$data[$i]}->{$i}+$num;
 }
 }

my @bp=("A","T","G","C");
	foreach my $key (@bp) {
		   print OUT3 "$key";   
			 for(my $i=0;$i<=29;$i++){
		     my  $total_r=$count{A}->{$i}+$count{C}->{$i}+$count{T}->{$i}+$count{G}->{$i};
			 if($total_r>0){
	   		 my    $rate=$count{$key}->{$i}/$total_r ;
		    printf OUT3 "\t%.3f",   "$rate";
		}
		} 
		print OUT3 "\n";

}

##############################################################################################3
  
#plot2 GC;
my $gc_pdf=$outdir."/GC.pdf";
my $gc_infile=$outdir."/GC_distribution";
my $r2=<<END;
gc<-read.table("$gc_infile",sep="\\t",header=TRUE)\n
 pdf("$gc_pdf")\n
 len<-length(gc[1,])\n
 gcc<-as.matrix(gc[,2:(len-1)])\n
 barplot(gcc, col=rainbow(4),main=c("Each site GC percent"), ylim=c(0,1), xlim=c(0,len+10), xlab="Clean read each site", ylab="%",font.lab=1)\n
 legend(len+3,0.8,legend=c("A","T","G","C"),fill=rainbow(4))\n
 dev.off()\n
END
  $R->run($r2);
  $R->stop();

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
