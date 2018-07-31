#!/usr/bin/perl -w
#Filename:
#Author: zhang hongmei
#Email: zhm.2009.happy@163.com
#Date: 2013/07/21;
#Modified: 
#Description: get miRNA expression from next generation sequencing data.
my $version=1.00;

#use strict;
use Getopt::Long;
use File::Basename;

if(@ARGV!=1){
	print "USAGE: perl \$0 config.properties\n";
}
my ($infile)=@ARGV;
my (%opt, $outdir, @sample);

open (OPT,"<$infile") or die "cannot open config.properties\n";
while(<OPT>){
	chomp;

	if(/RAW_DATA_PATH\s*=\s*(\S+)/)      {  $opt{'input'}=$1;}
	if(/RESULT_PATH\s*=\s*(\S+)/)		 {	$opt{'output'}=$1;}
	if(/FORMAT\s*=\s*(\S+)/)			 {	$opt{'format'}=$1;}
	if(/SAMPLE_LIST\s*=\s*(\S+)/)		 {  $opt{'sample'}=$1;}
	if(/GENOME_REFERENCE_FILE\s*=\s*(\S+)/) {  $opt{'genome_ref'}=$1;}
	if(/MASK_GENOME\s*=\s*(\S+)/) { $opt{'mask_genome'}=$1;}
	if(/NO_MIRNA_REFERENCE_FILE\s*=\s*(\S+)/) { $opt{'nomirna'}=$1;}
	if(/NO_MIRNA_REFERENCE_LENGTH_FILE\s*=\s*(\S+)/) { $opt{'nomirna_len'}=$1;}

	#miRNA info;
	if(/MAT_MIRNA\s*=\s*(\S+)/)            { $opt{'mat_txt'}=$1; }
	if(/PRE_MIRNA\s*=\s*(\S+)/)            { $opt{'pre_txt'}=$1; }
	
	#compare list
	if(/COMPARE_LIST\s*=\s*(.*)/)			{$opt{'compare_list'}=$1;}
	#path;
	if(/PERL_PATH\s*=\s*(\S+)/)            { $opt{'perl'}=$1;}
}
close OPT;

###put sample list into @sample
open(IN, "<$opt{'sample'}") || die "can not open $opt{'sample'}\n";
@sample = <IN>;
close IN;
$comp_dir=$opt{'output'}."/comparison/";
mkdir $comp_dir;

for(my $i=0; $i<=$#sample; $i++){
	chomp $sample[$i];
	my $sample_now=$sample[$i];
	$outdir[$i]=$opt{'output'}."/$sample[$i]"."/";
	print $outdir[$i],"\n";
	mkdir $outdir[$i];
	########step 0 adapter trim;
	my $outdir0=$outdir[$i]."data_filter/";
	mkdir $outdir0;
	my $outdir0_1=$outdir0."R1/";
	my $outdir0_2=$outdir0."R2/";
	mkdir $outdir0_1;
	mkdir $outdir0_2;
	`perl $opt{'perl'}"Adapter_trim_for_R1.pl" -i  $opt{'input'}$sample[$i]"_1.fq" -f $opt{'format'} -o $outdir0_1"/R1_clean.fa"`;
	`perl $opt{'perl'}"Adapter_trim_for_R2.pl" -i  $opt{'input'}$sample[$i]"_2.fq" -f $opt{'format'} -o $outdir0_2"/R2_clean.fa"`;
	`perl $opt{'perl'}"merge_same_reads_of_R1_R2.pl" -r1 $outdir0_1"/R1_clean.fa" -r2 $outdir0_2"/R2_clean.fa" -o $outdir0"/clean.fa"`;
	########step 1 genome mapping, 1 mismatch;
	my	$outdir1=$outdir[$i]."genome_mapping/";
	mkdir $outdir1;
	my	$samfile=$outdir1."sequence.sam";
	my	$saifile=$outdir1."sequence.sai";
	open OUT,">$samfile";
	open OUT2,">$saifile";
	close OUT;
	close OUT2;
	`perl $opt{'perl'}"genomapping_mis1.pl" -s $opt{'perl'} -sam $samfile -sai $saifile -genome $opt{'genome_ref'}  -queryseq $outdir0"/clean.fa"  -queue t `;
	
	########step 2 non_miRNA mapping, 1 mismatch;
	my $outdir2=$outdir[$i]."non-miRNA_mapping/";
	mkdir $outdir2;
	$opt{'bp'}=$outdir1."alignedTogenome.blastparsed";
	$opt{'fa'}=$outdir1."alignedTogenome.fa";
	`perl $opt{'perl'}"alignToNo_micRNA_zhm.pl" -s $opt{'perl'} -bp $opt{'bp'} -bf $opt{'fa'} -o $outdir2 -nomiRNA $opt{'nomirna'} -lf $opt{'nomirna_len'} -queue t`;

	########step 3 miRNA expression;
	my $outdir3=$outdir[$i]."miREXpress_edit/";
	mkdir $outdir3;
	open(IN, "<$outdir2"."mirna_candidate.fa");
	open(OUT,">$outdir2"."mirna_candidate.txt");
	while (my $aline=<IN>)
	{
		chomp $aline;
		if($aline=~/>.*_x(\d+)/)
		{
			print OUT "$1\t";
		}
		else
		{
			print OUT "$aline\n";
		}															
	}
	close IN;
	close OUT;
	##### 0 mismatch map to pre_miRNA
	`alignmentSIMD -r $opt{'pre_txt'} -i $outdir2"mirna_candidate.txt"  -u 10 -o $outdir3`; 
	`analysis -r  $opt{'pre_txt'}  -m $opt{'mat_txt'} -d $outdir3 -o sample_alignment_result -t sample_expression_profile`;
	##### 1 mismatch map nohit.fa to mask_genome
	my $infile2 = $outdir3."nohit";
	my $outfile2 = $outdir3."nohit.fa";
	
	$outdir4 = $outdir3."bowtie/";
	mkdir $outdir4;

	open(IN, "<$infile2") || die $!;
	open(OUT, ">$outfile2") || die $!;
	my $i = 0;
	my (@arr2,$num2,$seq2,$line2);
	while($line2 = <IN>){
		chomp $line2;
		@arr2 = split(/\t/, $line2);
		$num2 = $arr2[0];
		$seq2 = $arr2[1];
		print OUT ">nohit_".$i."_x".$num2."\n";
		print OUT $seq2."\n";
		$i++;
	}
	close IN;
	close OUT;

	`bowtie $opt{'mask_genome'} -f $outfile2 $outdir4"mapfile" --un $outdir4"unmapped" -a --best --strata -v 1`;
	
	#change unmapped nohit reads in to txt format which can be the input of mirexpress
	my $outdir5 = $outdir4."mirexpress_mis1/";
	mkdir $outdir5;

	my $infile3 = $outdir4."unmapped";
	my $outfile3 = $outdir4."unmapped.txt";

	open(IN, "<$infile3") || die $!;
	open(OUT, ">$outfile3") || die $!;
	my ($line3,$seq3,$num3);
	while($line3 = <IN>){
		$seq3 = <IN>;
		if($line3 =~ />nohit_\d+_x(\d+)\n/){
			$num3 = $1;
		}
		print OUT "$num3\t$seq3";
	}
	close IN;
	close OUT;
	
	#### 1_mismatch map the rest of nohit to pre_miRNA
	`alignmentSIMD -r $opt{'pre_txt'} -i $outfile3 -t 0.95  -u 10 -o $outdir5`; 
	`analysis -r  $opt{'pre_txt'}  -m $opt{'mat_txt'} -d $outdir5 -o sample_alignment_result -t sample_expression_profile`;
	#### plus the 0 mismatch and 1 mismatch result
	`cat $outdir3"/sample_alignment_result" $outdir5"/sample_alignment_result" >> $outdir3"tmp.alignment" `;
	`perl $opt{'perl'}"stat_miRNA_reads_by_alignment_result.pl" -i $outdir3"tmp.alignment" -o $outdir3"miRNA_reads"`;

	#### calculate miRNA RPM 
	$map_genome=$outdir1."map_to_genome.count";
	open(IN,"<$map_genome") or die $!;
	while(<IN>)
	{
		if(/^total\s+\d+\s+(\d+)\s+\d+/)
		{
			$total_num=$1;
		}
	}
	close IN;
	$rpm_out = $outdir3.$sample_now."_reads_rpm";
	print $sample_now;
	`perl $opt{'perl'}"calculate_miRNA_RPM.pl" -reads $outdir3"miRNA_reads" -seq_depth $total_num -rpm $rpm_out `;
	`cp $rpm_out $comp_dir `;
}

#### differential expression analysis, Poisson's distribution
my @compare_list=split(";",$opt{'compare_list'});
chdir $comp_dir;
#system("pwd");
foreach my $list (@compare_list)
{
	my @lists=split(",",$list);
	my $sample1=$lists[0]."_reads_rpm";
	my $sample2=$lists[1]."_reads_rpm";
	`perl $opt{'perl'}"ExpDiff.pl" -s1 $sample1 -s2 $sample2 -tpm 10 -outdir ./ -perl_path $opt{'perl'}`
}
