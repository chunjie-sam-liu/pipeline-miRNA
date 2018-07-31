#!/usr/bin/perl
#example : perl divide_conquer.pl xxx.fastq
	#yread from fastq file 
        $file =@ARGV[0];
        open(FILE,"<$file");
        @seq_list = <FILE>;
        for ($i=1;$i<=$#seq_list;$i+=4){
                push(@result,$seq_list[$i]);
		
        }
#	print "@result\n";

	#remove duplicates
        @result = grep { ++$hash{$_} < 2 } @result;
	#print "@result\n"; 


	#regexp to find longest common sequence
	sub regexp(){
		$subarray=shift ;
		@subarray=@$subarray;
		$can_total = join ("" , @subarray);
	       # my @substr = sort { length($b) <=> length($a) } $can_total =~ /(?=^.*$\n)*?.*?(.{8,})(?=.*\n)+.*\1/g;
	      my @substr = sort { length($b) <=> length($a) } $can_total=~ /(.{8,})(?=.*\n.*\1)/g;
	        my @res = grep { length == length $substr[0] } @substr;
		@res = grep { ++$hash{$_} < 2 } @res;
		return \@res;
	}
        sub find_adapter(){
		$sub_adapter=shift;
		$sub_adapter=$$sub_adapter;
	#	@seq_result=();
		foreach (@result){
			if ($position=index($_,$sub_adapter))
			{
				$single=substr($_,$position);
				push(@seq_result,$single);
#				print ;
			}
		}
		@seq_result2= sort {length ($b) <=> length ($a)} @seq_result;
		print $seq_result2[0];
	}

	#main function
	sub lcs(){
		$seq = shift;
		@seq=@$seq;	
	#	print "@$seq\n";
		for ($i=0;$i<$#seq;$i+=2){
			@subarray=@seq[$i..($i+1)];
		#	print "@subarray\n";
#		print @{&regexp(\@subarray)};	
			push(@tmp,@{&regexp(\@subarray)});
		}
		foreach (@tmp){
			$count{$_}++;
		}
		@seq_key=sort {$count{$b} <=> $count{$a}} keys %count;
		if ($count{$seq_key[0]}/($#result+1)>0.8)

		{
		#	print "$seq_key[0]\n";
			&find_adapter(\$seq_key[0]);
		}
		if($#tmp<=10)
		{
		#	print "$seq_key[0]\n";;
			&find_adapter(\$seq_key[0]);
		}
		else {
			@tmp = sort {$a <=> $b} @tmp;
			while (@tmp)
			{
				push (@next_circle,shift @tmp);
				push (@next_circle,pop @tmp);
			}
			&lcs(\@tmp)
		}
		@seq_key=();
	}
		
		
	&lcs(\@result);	
