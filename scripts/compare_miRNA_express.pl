#!/usr/bin/perl 
#Author: Gong Jing 2011-12-5; 
#Description: used to compare miRNA expression;
=infile microRNA.list format
#ID_1   tag     chrom   strand  pre_start       pre_end region  m_strat m_end   s_start s_end   mature  star    pre     id_2    m_exp   s_exp   m_exp_total     s_exp_total
hsa-miR1093     Conserve        chr10   +       120506472       120506523       1       22      27      52      AGGAGGTGAGCTGTGGGCCGGC  ATTACCACCTGAGCCCCACCTCCTGT      AGGAGGTGAGCTGTGGGCCGGCGAGCATTACCACCTGAGCCCCACCTCCTGT    smo-miR1093;tca-miR-3900-3p;hsa-miR-3689d;mtr-miR5254   19      0       19	0

usage: $0   result_file1 result_file2 -------
please save the result file in the same dirctory. first parameter is the dirctory, following parameters are result dirctory;

=cut
use File::Basename;

my @files=@ARGV; 
my $dir=dirname ($files[0]);
print "$dir\n";
$dir=~s/\/$//;
my $outdir=$dir."\/compare_result\/";                          # save the result in the /compare_result/ dirctory;
if( ! -e $outdir){
mkdir  $outdir ; 
}

my $i=0;
my (%data,%mirs);
my ($infile,$infile2,$mir,$tag,$chr,$mat,$star,$m_ex,$s_ex);
print "$#files\t$files[0]\t$files[1]\n";
for($i=0;$i<=$#files;$i++){
	$files[$i]=~s/\/$//;
	$infile=$files[$i]."\/target\/microRNA.list";
	$infile2=$files[$i]."\/filter\/sequence.stat";
	open(IN,"<$infile") or die "cant open $infile\n"; 
	while (<IN>){
		chomp;
		($mir,$tag,$chr,$strand,$pre_s,$pre_e,$mat,$star,$pre,$m_ex,$s_ex)=(split/\t/)[0,1,2,3,4,5,10,11,12,16,17];
		$miRNA1 =$mat;
		$miRNA2	=$star;
		$data[$i]{$miRNA1}=$m_ex;
		$data[$i]{$miRNA2}=$s_ex;
		$mirs{$miRNA1}="";
		$mirs{$miRNA2}="";
#	print "$data[$i]{$mir}{'chr'}\t $data[$i]{$mir}{'mat'}\t$data[$i]{$mir}{'star'}\n";
	}
	close IN;

	open(IN2,"<$infile2") or die "cannot open $infile2\n";
	while(<IN2>){
		if($_=~/polyA:\s*(\d+)/){$data[$i]{'total'}=$1;}           # count the total number of each data set;
	}
}

my $outfile=$outdir."miRDEEP_predicted_result\.xls"; 
open OUT,">$outfile" or die "cannot open $outfile\n"; 

#group the total miRNAs; eg:  a miRNA in  0,1,2,3,or 4 datasets;

foreach my  $mir2 (keys %mirs){
	
	my $j;
	my $n=0;
	for($j=0;$j<=$#files;$j++){
		if(exists  $data[$j]{$mir2}){
			$n++;
		}
	}
	$mir_group{$n}{$mir2}="";
}	

	for(my $n=$#files;$n>=0;$n--){ 
		foreach my $mir3 (keys %{$mir_group{$n}}){
			print OUT "$n\t$mir3\t";
		for(my $j=0;$j<=$#files;$j++){
		my	$read=$data[$j]{$mir3};
		my	$rpm_m=sprintf ("%8.2f",($read/$data[$j]{'total'}*1000000));
		#	my	$rpm_s=sprintf ("%8.2f",($read_s/$data[$j]{'total'}*1000000));
			print OUT "$read\t$rpm_m\t";
		}
		print OUT "\n";
	}
}



