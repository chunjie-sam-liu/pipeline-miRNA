#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $usage=
"Usage: $0 <in1.species.name> <in2.prediction> <in3.RNA.fasta> <in4.RNA.blastparsed><in5.GeneListChrom> <out1.conserved.xls> <out2.novel.xls> <out3.remained.fasta> <out4.loop.fasta><out5.novel_aln.result> 
 <out6.stat> <out7.stat>";

my $species = $ARGV[0] or die $usage;

open (INprediction, "< $ARGV[1]") or die $usage;
open (INfastaRNA, "< $ARGV[2]") or die $usage;
open (InblastparsedRNA, "< $ARGV[3]")or die $usage;
open (InGeneListChrom,"< $ARGV[4]")or die $usage;
open (OUTconserved, "> $ARGV[5]") or die $usage;
open (OUTnovel, "> $ARGV[6]") or die $usage;
open (OUTremained, "> $ARGV[7]") or die $usage;
open (OUTloop, "> $ARGV[8]") or die $usage;
open (OUTstat,"> $ARGV[9]")or die $usage;
open (OUTaln,">$ARGV[10]")or die $usage;
open (OUTmiRNASPE,">$ARGV[11]")or die $usage;
open (OUTNovelGff3,">$ARGV[12]") or die $usage ;

my ($a, @arr1, @arr2, $miRNA_ID,$miRNA_ID_checked,$miRNA,@mature,@pre, @pre_struct,@pos, @pri,@pri_struct,@star,$star_seq,$mature_seq,$pri_miRNA,$pre_miRNA, @rna, $count);

my ($b, $name, $seq, %hash, %mapped, %loop, @id, $num);
my ($c,@bp,%chrome,%start,%end,%strand);
my ($d,@genelist,%genelist_chrom);
my ($N_conserve, $N_conserve_mature, $N_conserve_mature_unique,$N_conserve_star,$N_conserve_star_unique);
my ($N_novel, $N_novel_mature, $N_novel_mature_unique,$N_novel_star,$N_novel_star_unique);
$N_conserve = $N_conserve_mature = $N_conserve_mature_unique = $N_conserve_star = $N_conserve_star_unique = 0;
$N_novel = $N_novel_mature = $N_novel_mature_unique = $N_novel_star = $N_novel_star_unique = 0;
my (@mature_arm,@mature_beg, @mature_end, @mature_query,@pri_beg,@pri_end,@star_beg, @star_end);
my (%rnaId,@conserved_miRNA,@conserved_miRNA_checked,@novel_miRNA,%miRNA_chrome,%miRNA_start,%miRNA_end,%miRNA_strand,%miRNA_region,%N_mature_total,%N_mature_unique,%N_star_total,%N_star_unique);
my $index;
my ($miRNA_chrom_start,$miRNA_chrom_end,$gene_chrom_start_temp,$gene_chrom_end_temp,$temp_gene_chrom);

$count=1;
#read remain-miRNA.fa
while ($b = <INfastaRNA>) {
	chomp $b;
	if ($b =~ /^>/) {
		$name = substr($b, 1);
		$b = <INfastaRNA>;
	    chomp $b;
	    $hash{$name} = $b;
	}	
}
close INfastaRNA;


#read remain-miRNA.blastparsed
while ($c=<InblastparsedRNA>) {
	chomp $c;
	@bp = split(/\t/, $c);
	$name = $bp[0];
	$chrome{$name}=$bp[3];
	@arr1=split(/\.\./, $bp[5]);
	$start{$name}=$arr1[0];
	$end{$name}=$arr1[1];
	$strand{$name}=$bp[9];
}
close InblastparsedRNA;

#read GeneListChrom
$index=1;
while ($d=<InGeneListChrom>) {
	chomp $d;
	
	@genelist=split(/\s/,$d);
	
	#strand
	$genelist_chrom{$index}[0]=$genelist[2];
	
	#start
	$genelist_chrom{$index}[1]=$genelist[3];
	
	#end
	$genelist_chrom{$index}[2]=$genelist[4];
	
	$index++;
}
close InGeneListChrom;


while ($a = <INprediction>) {

	until ($a =~ /^score_nucleus/) 
	{
		$a = <INprediction>;
		if (eof) {
		last;
	    }
	}
	if (eof) {
		last;
	}
	chomp $a;

	if ($a =~ /^score_nucleus	3/) {

		$N_conserve++;
		print OUTconserved "#miRNA_ID\t#pre_miRNA\t#pre_structure\t#mature_arm\t#mature_start\t#mature_end\t#star_start\t#star_end\n";

		$a = <INprediction>;
		chomp $a;
		@arr1 = split(/\t/, $a);
		
		$miRNA_ID_checked=$arr1[0];
		
		foreach  $miRNA(@arr1) {
			@arr2 = split(/\-/,$miRNA);
			
			if($arr2[0] eq $species){
				@arr2 = split(/\-/,$miRNA_ID_checked);
				$miRNA_ID_checked=~ s/$arr2[0]/$species/;
				last;
			}
		}
		
	
		$miRNA_ID = $arr1[0];
		@arr2 = split(/\-/, $miRNA_ID);
		$miRNA_ID =~ s/$arr2[0]/$species/;

		$miRNA = $arr1[0];
		$miRNA =~ s/$arr2[0]\-//;
		
		#first line:miRNA_id
		
		$conserved_miRNA_checked[$N_conserve]=$miRNA_ID_checked;
		if (!defined $rnaId{$miRNA_ID}) {
			$rnaId{$miRNA_ID} = 0;
			$conserved_miRNA[$N_conserve]=$miRNA_ID;
			print OUTmiRNASPE "$miRNA\t";
			print OUTconserved "$miRNA_ID\t";
			print OUTaln ">$miRNA_ID\n";
		}else{
			$rnaId{$miRNA_ID}++;
			$conserved_miRNA[$N_conserve]="$miRNA_ID-$rnaId{$miRNA_ID}";
			print OUTmiRNASPE "$miRNA-$rnaId{$miRNA_ID}\t";
			print OUTconserved "$miRNA_ID-$rnaId{$miRNA_ID}\t";
			print OUTaln ">$miRNA_ID-$rnaId{$miRNA_ID}\n";
		}
		foreach  $miRNA(@arr1) {
			@arr2 = split(/\-/,$miRNA);
			print OUTmiRNASPE "$arr2[0]";
			print OUTmiRNASPE ";";
		}
		print OUTmiRNASPE "\n";

		$N_mature_total{$conserved_miRNA[$N_conserve]}=$N_mature_unique{$conserved_miRNA[$N_conserve]}=$N_star_total{$conserved_miRNA[$N_conserve]}=$N_star_unique{$conserved_miRNA[$N_conserve]}=0;

		do{
			$a = <INprediction>;

		}until ($a=~ /^mature_arm/);
		chomp $a;
		@mature_arm = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@mature_beg = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@mature_end = split(/\t/,$a);		

		#mature_query
		$a = <INprediction>;
		chomp $a;
		@mature_query = split(/\t/,$a);
		$miRNA_chrome{$conserved_miRNA[$N_conserve]}=$chrome{$mature_query[1]};
		$miRNA_start{$conserved_miRNA[$N_conserve]}=$start{$mature_query[1]};
		$miRNA_end{$conserved_miRNA[$N_conserve]}=$end{$mature_query[1]};
		$miRNA_strand{$conserved_miRNA[$N_conserve]}=$strand{$mature_query[1]};

		#mature_seq
		do{
			$a=<INprediction>;
			if(eof){
				last;
			}
		}until ($a=~ /^mature_seq/);
		chomp $a;
		@mature = split(/\t/,$a);
		$mature_seq = $mature[1];
		
		do {
			$a = <INprediction>;
		} until ($a =~ /^pre_seq/) ;
		chomp $a;
		@pre = split(/\t/, $a);
		$pre_miRNA = $pre[1];
		print OUTconserved "$pre_miRNA\t";
		$a=<INprediction>;
		chomp $a;
		@pre_struct = split(/\t/,$a);
		print OUTconserved "$pre_struct[1]\t";
		print OUTconserved "$mature_arm[1]\t";
		print OUTconserved "$mature_beg[1]\t";
		print OUTconserved "$mature_end[1]\t";

		#pri_beg
		do{
			$a = <INprediction>;
		}until ($a =~ /^pri_beg/) ;
		chomp $a;
		@pri_beg = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@pri_end = split(/\t/,$a);
		
		#pri_seq
		do{
			$a=<INprediction>;
		}until ($a =~ /^pri_seq/);
		chomp $a;
		@pri = split(/\t/,$a);
		$pri_miRNA = $pri[1];
		print OUTaln "$pri_miRNA\t";
		if($rnaId{$miRNA_ID}==0){
			print OUTaln "$miRNA_ID\t";
		}
		else{
			print OUTaln "$miRNA_ID-$rnaId{$miRNA_ID}\t";
		}
		
		print OUTaln $pri_end[1]-$pri_beg[1]+1;
		print OUTaln "\n";

		#pri_struct
		$a = <INprediction>;
		chomp $a;
		@pri_struct = split(/\t/,$a);
		print OUTaln "$pri_struct[1]\n";
		
		#mature_seq
		$index = 1;
		do{
			print OUTaln "*";
			++$index;
		}until ($index == $mature_beg[1]);
		print OUTaln "$mature_seq";
		$index = $mature_end[1];
		do{
			print OUTaln "*";
			++$index;
		}until ($index == $pri_end[1]);


		print OUTaln "\tmature\t";
		print OUTaln $mature_end[1]-$mature_beg[1]+1;
		print OUTaln "\n";
				
		#star_beg
		do {
			$a = <INprediction>;
		} until ($a =~ /^star_beg/) ;
		chomp $a;
		@star_beg = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@star_end = split(/\t/,$a);
		
		print OUTconserved "$star_beg[1]\t";
		print OUTconserved "$star_end[1]\n";
		do {
			$a = <INprediction>;
		}until ($a =~ /^star_seq/);
		chomp $a;
		@star = split(/\t/,$a);
		$star_seq = $star[1];
		#star_seq
       		$index = 1;
		do{
			print OUTaln "*";
			++$index;
		}until ($index == $star_beg[1]);
		print OUTaln "$star_seq";
		$index = $star_end[1];
		do{
			print OUTaln "*";
			++$index;
		}until ($index == $pri_end[1]);

		print OUTaln "\tstar\t";
		print OUTaln $star_end[1]-$star_beg[1]+1;
		print OUTaln "\n";
				
		do {
			$a = <INprediction>;
		} until ($a =~ /^star_struct/) ;

		$a = <INprediction>;
		chomp $a;

		do {

		    @rna = split(/\t/, $a);
			@id = split(/\_/, $rna[0]);
			$num = substr($id[$#id],1);
			@pos = split(/\.\./, $rna[5]);

			if ($mature_end[1]<=$star_beg[1]) {#mature seq in first arm 
				if ($pos[1] <= $mature_end[1]) {					
				    
				#print out query sequences
					
					$index=1; 
					while($index<$pos[0]){
						print OUTaln "-";
						++$index;
					}
		
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";

					if (! defined $mapped{$rna[0]}) {
						$N_conserve_mature_unique ++;
						$N_conserve_mature += $num;	#conserve_mature total num 							
					}
					$N_mature_total{$conserved_miRNA[$N_conserve]}+=$num;
					$N_mature_unique{$conserved_miRNA[$N_conserve]}++;
			    }elsif($star_beg[1] <= $pos[0]){						
					$index=1;
					while($index < $pos[0]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";
					
					if (! defined $mapped{$rna[0]}) {
						$N_conserve_star_unique ++;
						$N_conserve_star += $num;	#conserve_mature unique num

					}
					$N_star_total{$conserved_miRNA[$N_conserve]}+=$num;
					$N_star_unique{$conserved_miRNA[$N_conserve]}++;
			    }else{
					if (! defined $loop{$rna[0]}) {
						$loop{$rna[0]} = 0;		
					}
				}
			}else{#mature seq in second arm
				if ($mature_beg[1] <= $pos[0]) {					
					#print out query sequences
					
					$index=1; 
					while($index <$pos[0]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";
					if (! defined $mapped{$rna[0]}) {
						$N_conserve_mature_unique ++;
						$N_conserve_mature += $num;	

					}
					$N_mature_total{$conserved_miRNA[$N_conserve]}+=$num;
					$N_mature_unique{$conserved_miRNA[$N_conserve]}++;
			    }elsif($pos[1] <= $star_end[1]){					
				    
					#print out query sequences
					
					$index=1; 
					while($index < $pos[0]){
						print OUTaln "-";
						++$index;
					}
				
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";
					if (! defined $mapped{$rna[0]}) {
						$N_conserve_star_unique ++;
						$N_conserve_star += $num;

					}
					$N_star_total{$conserved_miRNA[$N_conserve]}+=$num;
					$N_star_unique{$conserved_miRNA[$N_conserve]}++;
			    }else{
					if (! defined $loop{$rna[0]}) {
						$loop{$rna[0]} = 0;		
					}
				}
			}			

			if (! defined $mapped{$rna[0]}) {
				$mapped{$rna[0]} = 0;			
			}
			else{
				$mapped{$rna[0]}++;
			}

			if (eof) {
				last;
			}

			$a = <INprediction>;
		    chomp $a;
			until ($a !~ /Minus/) {
				$a = <INprediction>;
		        chomp $a;
			}

		} until (length($a)<1) ;

	}else{

		$N_novel++;
		print OUTnovel "#miRNA_ID\t#pre_miRNA\t#pre_structure\t#mature_arm\t#mature_start\t#mature_end\t#star_start\t#star_end\n";
		$novel_miRNA[$N_novel]="$species-miR-c-$count";
		print OUTnovel "$species-miR-c-$count\t";

		#print align message to novel_aln.result
		#first line:miRNA_id
		print OUTaln ">$species-miR-c-$count\n";
	
		$N_mature_total{$novel_miRNA[$N_novel]}=$N_mature_unique{$novel_miRNA[$N_novel]}=$N_star_total{$novel_miRNA[$N_novel]}=$N_star_unique{$novel_miRNA[$N_novel]}=0;

		do{
			$a = <INprediction>;

		}until ($a=~ /^mature_arm/);
		chomp $a;
		@mature_arm = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@mature_beg = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@mature_end = split(/\t/,$a);

		#mature_query
		$a = <INprediction>;
		chomp $a;
		@mature_query = split(/\t/,$a);
		$miRNA_chrome{$novel_miRNA[$N_novel]}=$chrome{$mature_query[1]};
		$miRNA_start{$novel_miRNA[$N_novel]}=$start{$mature_query[1]};
		$miRNA_end{$novel_miRNA[$N_novel]}=$end{$mature_query[1]};
		$miRNA_strand{$novel_miRNA[$N_novel]}=$strand{$mature_query[1]};

		$miRNA_region{$novel_miRNA[$N_novel]}="Intergenic";
		foreach my $iindex (keys %genelist_chrom) {
			#print $genelist_chrom{$iindex}[0];
			if($genelist_chrom{$iindex}[0] eq $strand{$mature_query[1]}&&$genelist_chrom{$iindex}[1]<=$start{$mature_query[1]} && $genelist_chrom{$iindex}[2]>=$end{$mature_query[1]}){
				$miRNA_region{$novel_miRNA[$N_novel]}="Intron";
				last;
			}
		}

        #mature_seq
		do{
			$a=<INprediction>;
			if(eof){
				last;
			}
		}until ($a=~ /^mature_seq/);
		chomp $a;
		@mature = split(/\t/,$a);
		$mature_seq = $mature[1];
		
		do {
			$a = <INprediction>;
			if (eof) {
				last;
			}
		} until ($a =~ /^pre_seq/) ;
		chomp $a;
		@pre = split(/\t/, $a);
		$pre_miRNA = $pre[1];
		print OUTnovel "$pre_miRNA\t";
		$a=<INprediction>;
		chomp $a;
		@pre_struct = split(/\t/,$a);
		print OUTnovel "$pre_struct[1]\t";
		print OUTnovel "$mature_arm[1]\t";
		print OUTnovel "$mature_beg[1]\t";
		print OUTnovel "$mature_end[1]\t";
		

		#pri_beg
		do{
			$a = <INprediction>;
		}until ($a =~ /^pri_beg/) ;
		chomp $a;
		@pri_beg = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@pri_end = split(/\t/,$a);
		
		#pri_seq
		do{
			$a=<INprediction>;
		}until ($a =~ /^pri_seq/);
		chomp $a;
		@pri = split(/\t/,$a);
		$pri_miRNA = $pri[1];
		print OUTaln "$pri_miRNA\t";
		print OUTaln "$species-miR-c-$count\t";
		print OUTaln $pri_end[1]-$pri_beg[1]+1;
		print OUTaln "\n";

		#pri_struct
		$a = <INprediction>;
		chomp $a;
		@pri_struct = split(/\t/,$a);
		print OUTaln "$pri_struct[1]\n";
		
		#mature_seq
		$index = 1;
		do{
			print OUTaln "*";
			++$index;
		}until ($index == $mature_beg[1]);
		print OUTaln "$mature_seq";
		$index = $mature_end[1];
		do{
			print OUTaln "*";
			++$index;
		}until ($index == $pri_end[1]);
		print OUTaln "\tmature\t";
		print OUTaln $mature_end[1]-$mature_beg[1]+1;
		print OUTaln "\n";
				
		#star_beg
		do {
			$a = <INprediction>;
		} until ($a =~ /^star_beg/) ;
		chomp $a;
		@star_beg = split(/\t/,$a);
		$a = <INprediction>;
		chomp $a;
		@star_end = split(/\t/,$a);

		print OUTnovel "$star_beg[1]\t";
		print OUTnovel "$star_end[1]\n";
		
		do {
			$a = <INprediction>;
		}until ($a =~ /^star_seq/);
		chomp $a;
		@star = split(/\t/,$a);
		$star_seq = $star[1];
		#star_seq
        $index = 1;
		while($index < $star_beg[1]){
			print OUTaln "*";
			++$index;
		}

		print OUTaln "$star_seq";
		$index = $star_end[1];
		while($index < $pri_end[1]){
			print OUTaln "*";
			++$index;
		}
		
		print OUTaln "\tstar\t";
		print OUTaln $star_end[1]-$star_beg[1]+1;
		print OUTaln "\n";
		
		do {
			$a = <INprediction>;
		} until ($a =~ /^star_struct/) ;
			
		$a = <INprediction>;
		chomp $a;
	
		do {

		    @rna = split(/\t/, $a);
			@id = split(/\_/, $rna[0]);
			$num = substr($id[$#id],1);
			@pos = split(/\.\./, $rna[5]);

			if ($mature_end[1]<=$star_beg[1]) {
				if ($pos[1] <= $mature_end[1]) {					
				    	
					#print out query sequences
					
					$index=1;
					while($index < $pos[0]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";
					
					if (! defined $mapped{$rna[0]}) {
						$N_novel_mature_unique++;
						$N_novel_mature += $num;	
					}
					$N_mature_total{$novel_miRNA[$N_novel]}+=$num;
					$N_mature_unique{$novel_miRNA[$N_novel]}++;
			    }elsif($star_beg[1] <= $pos[0]){					
				    
					#print out query sequences
					
					$index=1;
					while($index < $pos[0]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";
					if (! defined $mapped{$rna[0]}) {
						$N_novel_star_unique ++;
						$N_novel_star += $num;	

					}
					$N_star_total{$novel_miRNA[$N_novel]}+=$num;
					$N_star_unique{$novel_miRNA[$N_novel]}++;
			    }else{
					if (! defined $loop{$rna[0]}) {
						$loop{$rna[0]} = 0;		
					}
				}
			}else{
				if ($mature_beg[1] <= $pos[0]) {					
					#print out query sequences
					
					$index=1; 
					while($index < $pos[0]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";
					if (! defined $mapped{$rna[0]}) {
						$N_novel_mature_unique++;
						$N_novel_mature += $num;

					}
					$N_mature_total{$novel_miRNA[$N_novel]}+=$num;
					$N_mature_unique{$novel_miRNA[$N_novel]}++;
			    }elsif($pos[1] <= $star_end[1]){					
				    	
					#print out query sequences
					
					$index=1;

					while($index < $pos[0]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "$hash{$rna[0]}";
					$index=$pos[1];
					while($index < $pri_end[1]){
						print OUTaln "-";
						++$index;
					}
					
					print OUTaln "\t";
					print OUTaln "$rna[0]\t";
					print OUTaln $pos[1]-$pos[0]+1;
					print OUTaln "\n";
					if (! defined $mapped{$rna[0]}) {
						$N_novel_star_unique++;
						$N_novel_star += $num;	

					}
					$N_star_total{$novel_miRNA[$N_novel]}+=$num;
					$N_star_unique{$novel_miRNA[$N_novel]}++;
			    }else{
					if (! defined $loop{$rna[0]}) {
						$loop{$rna[0]} = 0;		
					}
				}
			}

			if (! defined $mapped{$rna[0]}) {
				$mapped{$rna[0]} = 0;			
			}
			else{
				$mapped{$rna[0]}++;
			}
			
			if (eof) {
				last;
			}

			$a = <INprediction>;
		    chomp $a;
			until ($a !~ /Minus/) {
				$a = <INprediction>;
		        chomp $a;
			}

		} until (length($a)<1) ;

		$count++;
	}
}


foreach my $key (keys %hash) {
	if (!defined $mapped{$key}) {
		print OUTremained ">$key\n";
		print OUTremained "$hash{$key}\n";
	}
}


foreach my $key (keys %hash) {
	if (defined $loop{$key}) {
		print OUTloop ">$key\n";
		print OUTloop "$hash{$key}\n";
	}
}

#print OUTstat

print OUTstat "conserved statics:\n";
print OUTstat "#miRNA_id\t#unique reads\t#total reads\t#chromosome\t#chrom_start\t#chrom_end\t#strand\t#miRNA_id_checked\n";
for($index=1;$index<=$N_conserve;$index++){
	print OUTstat "$conserved_miRNA[$index]\t";
	print OUTstat $N_mature_unique{$conserved_miRNA[$index]}+$N_star_unique{$conserved_miRNA[$index]};
	print OUTstat "\t";
	print OUTstat $N_mature_total{$conserved_miRNA[$index]}+$N_star_total{$conserved_miRNA[$index]};
	print OUTstat "\t";
	print OUTstat "$miRNA_chrome{$conserved_miRNA[$index]}\t";
	print OUTstat "$miRNA_start{$conserved_miRNA[$index]}\t";
	print OUTstat "$miRNA_end{$conserved_miRNA[$index]}\t";
	print OUTstat "$miRNA_strand{$conserved_miRNA[$index]}\t";
	print OUTstat "$conserved_miRNA_checked[$index]\n";
		
}

print OUTstat "\nnovel statics:\n";
print OUTstat "#miRNA_id\t#unique reads\t#total reads\t#chromosome\t#chrom_start\t#chrom_end\t#strand\t#chrom_region\n";
for($index=1;$index<=$N_novel;$index++){
	print OUTstat "$novel_miRNA[$index]\t";
	print OUTstat $N_mature_unique{$novel_miRNA[$index]}+$N_star_unique{$novel_miRNA[$index]};
	#print OUTstat $N_star_unique{$novel_miRNA[$index]};
	print OUTstat "\t";
	print OUTstat $N_mature_total{$novel_miRNA[$index]}+$N_star_total{$novel_miRNA[$index]};
	print OUTstat "\t";
	print OUTstat "$miRNA_chrome{$novel_miRNA[$index]}\t";
	print OUTstat "$miRNA_start{$novel_miRNA[$index]}\t";
	print OUTstat "$miRNA_end{$novel_miRNA[$index]}\t";
	print OUTstat "$miRNA_strand{$novel_miRNA[$index]}\t";
	print OUTstat "$miRNA_region{$novel_miRNA[$index]}\n";
}
print OUTstat "\ngeneral statics:\n";
print OUTstat "#num\t#miRNA\t#mature_unique\t#mature_total\t#star_unique\t#star_total\n";
print OUTstat "#conserved_miRNA\t$N_conserve\t$N_conserve_mature_unique\t$N_conserve_mature\t$N_conserve_star_unique\t$N_conserve_star\n";
print OUTstat "#novel_miRNA\t$N_novel\t$N_novel_mature_unique\t$N_novel_mature\t$N_novel_star_unique\t$N_novel_star\n";


##########print gff3
print OUTNovelGff3 "##gff-version 3\n";
for($index=1;$index<=$N_novel;$index++){
	print OUTNovelGff3 "$miRNA_chrome{$novel_miRNA[$index]}\t";
	print OUTNovelGff3 "novel\tmRNA\t";
	print OUTNovelGff3 "$miRNA_start{$novel_miRNA[$index]}\t$miRNA_end{$novel_miRNA[$index]}\t.\t$miRNA_strand{$novel_miRNA[$index]}\t.\tName=$novel_miRNA[$index];\n";
}

for($index=1;$index<=$N_conserve;$index++){
	print OUTNovelGff3 "$miRNA_chrome{$conserved_miRNA[$index]}\t";
	print OUTNovelGff3 "conserve\tmRNA\t";
	print OUTNovelGff3 "$miRNA_start{$conserved_miRNA[$index]}\t$miRNA_end{$conserved_miRNA[$index]}\t.\t$miRNA_strand{$conserved_miRNA[$index]}\t.\tName=$conserved_miRNA[$index];\n";
}


close INprediction;
close OUTconserved;
close OUTnovel;
close OUTremained;
close OUTloop;
close OUTstat;
close OUTaln;
close OUTmiRNASPE;
close OUTNovelGff3;
