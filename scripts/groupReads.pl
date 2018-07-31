#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2010-01
#Modified:
#Description:  
my $version=1.00;

use strict;
use Getopt::Long;
use File::Basename;

my %opts;
GetOptions(\%opts,"input=s","o=s","minNr:i","tag=s","qual:s","h");
if (!(defined $opts{input} and defined $opts{o} and defined $opts{tag} ) || defined $opts{h}) { #necessary arguments
&usage;
}
my ($qualT,$qualV);
if (defined $opts{'qual'}) {
	my @temp=split /:/,$opts{'qual'};
	$qualT=$temp[0];
	$qualV=$temp[1];
}
else{
$qualT="mean";
$qualV=0;
}
my @groups=split/\#/,$opts{'input'};
my $minNr=defined $opts{'minNr'}? $opts{'minNr'}:0;
my $outFile=$opts{'o'}.".txt";
my $tag=$opts{'tag'};
my @outputData;
my $fileout;
if($#groups >= 1){

	my @files1 = split(/,/,$groups[0]);
	my @files2 = split(/,/,$groups[1]);
	my (%count1, %count2, %allReads,$c1,$c2);
	$c1 = &getCountSum(\@files1,\%count1,$tag);
	$c2 = &getCountSum(\@files2,\%count2,$tag);
	&getAllReads(\%allReads, \%count1, \%count2); # get a "dummy" hash which holds all reads
	&getLog2(\%allReads, \%count1, \%count2, $c1, $c2);
	
	while(my @a = each(%allReads)){
		if($count1{$a[0]} and $count2{$a[0]}){
			my $sum = ($count1{$a[0]} + $count2{$a[0]})/2;
			if($sum >= $minNr){
				push @outputData,"$a[0],$sum,$a[1]";		
			}
		}
		elsif($count1{$a[0]}){
			my $sum = $count1{$a[0]}/2;
			if($sum >= $minNr){
				push @outputData,"$a[0],$sum,$a[1]";		
			}						
		}
		elsif($count2{$a[0]}){
			my $sum = $count2{$a[0]}/2;
			if($sum >= $minNr){
				push @outputData,"$a[0],$sum,$a[1]";		
			}
		}
	}
	  @outputData = &Sort_Table(cols=>3,field=>2,sorting=>"descending",
	  structure=>"csv",separator=>"\,",data=>\@outputData); 
	
	&writeOutput($outFile,\@outputData);
	my @final;
	push @final,@files1;
	push @final,@files2;
	
	if ($tag eq "solid") {$fileout = $opts{'o'}.".fa";}
	else{$fileout=$opts{'o'}.".fa";}
	&getname(\@final,$fileout,\%allReads,$tag);
}
# just one group --> "normal mode"
elsif($#groups == 0){

	my @files = split(/,/,$groups[0]);
	my %count;
	my $count = &getCountSum(\@files,\%count,$tag);
	
	while(my @a = each(%count)){
		if($a[1] >= $minNr){
			push @outputData,"$a[0],$a[1]";			
		}
	}
	
	@outputData = &Sort_Table(cols=>2,field=>2,sorting=>"descending",
	 structure=>"csv",separator=>"\,",data=>\@outputData); 
	&writeOutput($outFile,\@outputData);
	
	if ($tag eq "solid") {$fileout = $opts{'o'}.".fa";}
	elsif($tag eq "solexa"){$fileout = $opts{'o'}.".fa";}
	else{die "Error!pless check -tag !\n";}
	&getname(\@files,$fileout,\%count,$tag);
}
else{
	print ("an error ocurred - Please see the command line format \n\n");
	die;

}

######################################################################
#subprogram
######################################################################
sub getname(){
	#`touch "$_[1]"`;
	#`rm "$_[1]"`;
	if ($_[3] eq "solid") {
		my @fi = @{$_[0]};
		for(my $i = 0; $i <= $#fi; $i++){
			#print "go for $fi[$i] \n";
			&getfasta($fi[$i],$_[1],$_[2]);
		}

		
	}
	elsif($_[3] eq "solexa"){
		my @fi = @{$_[0]};
		for(my $i = 0; $i <= $#fi; $i++){
			#print "go for $fi[$i] \n";
			&getfastq($fi[$i],$_[1],$_[2]);
		}
	}
	else{
		print "an error ocurred - Please check the parameter:tag \n\n";
		die;
	}
}
sub getfasta(){
	open IN,"<$_[0]" or die "could not open $_[0]!\n";
	open OUT ,">>$_[1]";
	my $name;my $ref;
	while (my $aline=<IN>) {
		if ($aline=~/^\#/) {next;}
		my @name=split/\s+/," $aline";
		$name=$name[1];
		$aline=<IN>;
		my @ref=split/\s+/," $aline";
		$ref=$ref[1];
		if (defined $_[2]->{$ref}) {
			print OUT $name,"\n",$ref,"\n";
			delete ($_[2]->{$ref});
		}
	}
	close OUT;	
}

sub getfastq(){
	open IN,"<$_[0]" or die "could not open $_[0]!\n";
	open OUT ,">>$_[1]";
	my ($name,$ref,$qv);
	while (my $aline=<IN>) {
		chomp $aline;
		my @name=split/:/,"$aline";
		$name=">";
		for (my $i=1;$i<@name;$i++) {$name .=$name[$i]."_";}
		$aline=<IN>;
		my @ref=split/\s+/," $aline";
		$ref=$ref[1];
		$aline=<IN>;
		$aline=<IN>;
		my @qv=split/\s+/," $aline";
		$qv=$qv[1];
		if (defined $_[2]->{$ref}) {
#			print OUT $name,"\n",$ref,"\n","+\n",$qv,"\n";
			print OUT $name,"x",$_[2]->{$ref},"\n",$ref,"\n";
			delete ($_[2]->{$ref});
		}
	}
	close OUT;	
#	close FA;
}

sub writeOutput(){
	
	open (O,">$_[0]") or die "could not open $_[0]";
	my $last = @{$_[1]};
	for(my $i = 0; $i < $last; $i++){
		my $str = $_[1]->[$i];
		$str =~ s/,/\t/g;
		print O $str."\n";
	}
	close(O);
}

sub getLog2{
	
	while(my @a = each(%{$_[0]})){
		
		if($_[1]->{$a[0]} and $_[2]->{$a[0]}){
			
			my $ratio = ($_[1]->{$a[0]}/$_[3])/($_[2]->{$a[0]}/$_[4]);
			$_[0]->{$a[0]} = log($ratio)/log(2);
		}
		# expression just in "cases"
		elsif($_[1]->{$a[0]}){
			$_[0]->{$a[0]} = 10;			
		}
		# expression just in "control"
		elsif($_[2]->{$a[0]}){
			$_[0]->{$a[0]} = -10;						
		}
		
	}
	
}


sub getAllReads(){
	
	foreach my $k (keys %{$_[1]}) {
		$_[0]->{$k}++;
	}

	foreach my $k (keys %{$_[2]}) {
		$_[0]->{$k}++;		
	}
	
}

sub getCountSum(){
	
	my $totCount = 0 ;
	my @fi = @{$_[0]};
	for(my $i = 0; $i <= $#fi; $i++){
		print "go to $fi[$i] \n";
		my $c = getCounts($fi[$i],$_[1],$_[2]);
		$totCount += $c;
	}
	return $totCount;
}

sub getCounts(){
		if ($qualT eq "mean" and $_[2] eq "solid") {
			return &hashReadsSolid($_[0],$_[1]);
		}
		elsif($qualT eq "mean" and $_[2] eq "solexa"){
			return &hashReadsMean($_[0],$_[1]);
		}
		elsif($qualT eq "min" and $_[2] eq "solexa"){
			return &hashReadsMin($_[0],$_[1]);			
		}
		else{
			print "something wrong with the quality input option \n\n\n";
			die;
			
		}
}
sub hashReadsSolid(){
	open(I,$_[0]) or die "could not open $_[0]";
	my $counter=0;
	while (my $z=<I>) {
		if ($z=~/^\#/) {next;}
		$z=<I>;
		my @f = split(/\s+/," $z");
		my $seq = $f[1];
		$counter++;
		$_[1]->{$seq}++;
	}
	close I;
	return $counter;
}

sub hashReadsMean(){
	
	my $counter = 0;
	open(I,$_[0]) or die "could not open $_[0]";
	while(my $z = <I>){
		
		
		$z = <I>;
		my @f = split(/\s+/," $z");
		my $seq = $f[1];
		my $header = <I>; # the quality header

		$z = <I>;
		@f = split(/\s+/," $z");
		my $qual = $f[1];
		my $q = &getMeanQuality($qual);
		if($q > $qualV){
			$counter++;
			$_[1]->{$seq}++;
		}
	}
	return $counter;
}


sub hashReadsMin(){
	
	my $counter = 0;
	open(I,$_[0]) or die "could not open $_[0]";
	while(my $z = <I>){
		
		$z = <I>;
		my @f = split(/\s+/," $z");
		my $seq = $f[1];
		my $header = <I>; # the quality header

		$z = <I>;
		@f = split(/\s+/," $z");
		my $qual = $f[1];

		my $q = &getMinQuality($qual);
		if($q >= $qualV){
			$counter++;
			$_[1]->{$seq}++;
		}
	}
	return $counter;
}

#
# This function gives back the mean Q-value of the sequence/quality string
sub getMeanQuality(){
	
	my @bases = split(//,$_[0]);
	my $sum = 0;
	for(my $i = 0; $i <= $#bases; $i++){
		my $num = ord($bases[$i]) - 64;
		$sum += $num;
	}
	
	return $sum/($#bases+1);
	
}

###
### This function gives back the Q-value of the worst base
sub getMinQuality(){
	
	my @bases = split(//,$_[0]);
	my $worst = 1000;
	for(my $i = 0; $i <= $#bases; $i++){
#		printf ("base: $bases[$i]  --> %d\n",ord($bases[$i]));
		my $num = ord($bases[$i]) - 64;
		if($num < $worst){
			$worst = $num;
		}
	}
	return $worst;
}
sub Sort_Table {
	# Get the args and put them into a Hash.
   	my (%arg) = @_;
	my $error = 0;


	# Subtract 1 for better readable Arrayfields ->
	# beginning count at 1 (not 0). ;)
	$arg{cols}--;
	$arg{field}--;

	if ($arg{structure} eq 'single') {
		# Array is not semicolon-separated and we must
		# convert it to semicolon-separated.
		@_ = ();
		my $i=0;
		while (defined ${$arg{data}}[$i] ne '') {
			my $tmp='';
			for (0..$arg{cols}) {
				$tmp .= "${$arg{data}}[$i+$_]";
				if ($_ != $arg{cols}) {
					$tmp .= "$arg{separator}";
				}
			}
			push(@_, $tmp);
			$i += $arg{cols} + 1;
		}
		@{$arg{data}} = @_;
	}

	my $use_warn = 0;
	# Turn warnings off, because we do first a '<=>' and if that
	# fails, we do a 'cmp' and then a warning comes up.
	# After sorting, we turn $^W to the same as before.
	if ($^W) {
		$use_warn = $^W;
		$^W = 0;
	}
	if ($arg{sorting} eq 'ascending') {
		# Sorting content ascending order.
		@{$arg{data}} =
			map { $_->[0] }
			sort {
				$a->[1] <=> $b->[1]
					||
				$a->[1] cmp $b->[1]
			}
			map { [ $_, (split(/$arg{separator}/))[$arg{field}] ] }
		@{$arg{data}};
	}
	elsif ($arg{sorting} eq 'descending') {
		# Sorting content descending order.
		@{$arg{data}} =
			map { $_->[0] }
			sort {
				$b->[1] <=> $a->[1]
					||
				$b->[1] cmp $a->[1]
			}
			map { [ $_, (split(/$arg{separator}/))[$arg{field}] ] }
		@{$arg{data}};
	}

	# Turn warnings to the same as before.
	if ($use_warn) {
		$^W = $use_warn;
	}

	# Return the sorted Array in the
	# same format as input.
	if ($arg{structure} eq 'csv') {
		return @{$arg{data}};
	}
	elsif ($arg{structure} eq 'single') {
		@_ = ();
		foreach (@{$arg{data}}) {
			push(@_, split(/$arg{separator}/));
		}
		return @_;
	}
} 
	




sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -input -o -tag -minNr -qual
options:

-input #input file 
	eg:  (file1_group1,file2_group1#file1_group2,file2_group2)      
	  Input Files: The perl script needs the s_L_sequence.txt       
	  file from the Gerald folder, being L the lane.                 
	  the different files of the same expriment are coma separated   
	  in case you have two groups (case/control) you can provide     
	  the files of the two groups separated by and #                 

-o  #output file prefix

-tag #[solid/solexa]

-minNr #integer
	default 0;
	the minimum number of copies a sequence read must have to be written out

-qual #reads filter
   eg:(min:value/mean:value)
   This parameter just for solexa reads.
   If the input files are solid and needs filter,please do filter first .
   
-h help
USAGE
exit(1);
}

