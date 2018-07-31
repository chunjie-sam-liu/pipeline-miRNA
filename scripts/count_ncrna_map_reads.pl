#!/usr/bin/perl -w
#Filename:
#Author: 
#Email: 
#Date: 
#Modified:
#Description:

use strict;
use warnings;

open INgenome, "< $ARGV[0]";
open INncrna, "< $ARGV[1]";
open INmirna, "< $ARGV[2]";
open OUT, "> $ARGV[3]";

my ($i,$a, @arr1, @arr2,$temp, $line, @arr,@counter);
for ($i=0;$i<6;$i++){
	$counter[$i]=0;
}

#genome mapping unique reads number
#genome mapping total reads number
my %hash_genome;
while($a = <INgenome>) 
{
	chomp $a;
	@arr1 = split/\t+/,$a;
	$hash_genome{$arr1[0]}=1;
}
close INgenome;

foreach my $key (keys %hash_genome)
{
	@arr2 = split/_x/,$key;
	$temp = $arr2[1];
	$counter[0]++;
	$counter[1]=$counter[1]+$temp;	
}

#ncrna mapping unique reads number
#ncrna mapping total reads number
my %hash_ncrna;
while($a=<INncrna>)
{
	chomp $a;
	if($a=~/^>/)
	{
		next;
	}
	else
	{
		$hash_ncrna{$a}=1;
	}
}	
close INncrna;

foreach my $key (keys %hash_ncrna)
{
		$counter[2]++;
		@arr=split/_x/,$key;
		$counter[3]=$counter[3]+$arr[1];
}


#mirna mapping unique reads number
#mirna mapping total reads number
my %hash_mirna;
while($a = <INmirna>) 
{
	chomp $a;
	@arr1 = split/\t+/,$a;
	$hash_mirna{$arr1[0]}=1;
}
close INmirna;

foreach my $key (keys %hash_mirna)
{
	@arr2 = split/_x/,$key;
	$temp = $arr2[1];
	$counter[4]++;
	$counter[5]=$counter[5]+$temp;
	
}


print OUT "unique genome reads: $counter[0]\n";
print OUT "total genome reads: $counter[1]\n";
print OUT "unique ncrna reads: $counter[2]\n";
print OUT "total ncrna reads: $counter[3]\n";
print OUT "unique mirna reads: $counter[4]\n";
print OUT "total mirna reads: $counter[5] \n";

close INgenome;
close INncrna;
close INmirna;
close OUT;
