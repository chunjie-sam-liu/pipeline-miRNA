#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:
#Description: É¾³ýmatched reads 
my $version=1.00;

use strict;
use Getopt::Long;
#use Statistics::R;

my %opts;
GetOptions(\%opts,"i=s","o=s","h");
if (!(defined $opts{i} and defined $opts{o} ) || defined $opts{h}) { #necessary arguments
&usage;
}

open IN,"<$opts{i}";
open OUT ,">$opts{o}";
my ($uniq,$total,$uniq_map,$total_map,$uniq_un,$total_un);
$uniq=$total=$uniq_map=$total_map=$uniq_un=$total_un=0;
while (my $aline=<IN>) {
	chomp $aline;
	next if($aline=~/^\@/);
	my @temp=split/\t/,$aline;
	$temp[0]=~/_x(\d+)$/;
	$uniq ++;
	$total +=$1;
	if ($temp[1]==0 || $temp[1]==16) {
		$total_map +=$1;
		$uniq_map++;
		next;
	}
	$uniq_un ++;
	$total_un +=$1;
}

print OUT "class\ttotal\tmapped\tunmapped\n";
print OUT "total\t$total\t$total_map\t$total_un\n";
print OUT "unique\t$uniq\t$uniq_map\t$uniq_un\n";
close OUT;

####plot ;
#my $R = Statistics::R->new() ;  
#  $R->startR ;
#  my $r=<<END;
#  outfile<-paste("$opts{o}","_total.pdf",sep="")
#    data<-read.table("$opts{o}",sep="\\t",header=TRUE)\n
#	  pdf(outfile)\n
#	  b<-data[1,3:4]
#	  c<-as.matrix(data[1,3:4])
#	  color<-c("orange","green")
#	  names(c)<-c(data[1,3],data[1,4])
#	  xrange<-c(-2,2)
#	   m<-paste("Mapped to genome:",round(data[1,3]/data[1,2]*100,2),"%")
#	    u<-paste("Unmapped to genome:",100-round(data[1,3]/data[1,2]*100,2),"%")
#		pie(c,col=color,xlim=xrange,main=c("Total clean reads map to genome"))
#		legend("bottom",legend=c(m,u),fill=color)
	 
#	  dev.off()\n

#	   outfile2<-paste("$opts{o}","_uniq.pdf",sep="")
#	          data<-read.table("$opts{o}",sep="\\t",header=TRUE)\n
#                       pdf(outfile2)\n
#		              b<-data[2,3:4]
#		               c<-as.matrix(data[2,3:4])
#		                      color<-c("orange","green")
#	                       names(c)<-c(data[1,3],data[1,4])
#             xrange<-c(-2,2)
#                m<-paste("Mapped to genome:",round(data[2,3]/data[2,2]*100,2),"%")
#	                      u<-paste("Unmapped to genome:",100-round(data[2,3]/data[2,2]*100,2),"%")
#                    pie(c,col=color,xlim=xrange,main=c("Uniq clean reads map to genome"))
#		                  legend("bottom",legend=c(m,u),fill=color)
#             dev.off()\n	
#END

#	   $R->run($r);
#	   $R->stop();
sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i -o
options:
-r sam file
-o output file
-h help
USAGE
exit(1);
}

