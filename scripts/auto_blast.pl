#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $usage=
"$0 file_query file_subject temp_out_directory perl_script_directory

This is a wrapper for NCBI blastn.
It takes as input a fasta query file and a fasta subject (database) file.
It outputs the alignments in blast_parsed format.

Options:
-a Only print out perfect alignments (full coverage, 100% identity)
-b Only print out perfect alignments and sort the lines to signature input to miRDeep
-c Only print out alignments with this or lower e-value (default 1).
";


#Input files
my $file_query=shift or die $usage;
my $file_subject=shift or die $usage;
my $tmpdir=shift or die $usage;
my $perldir=shift or die $usage;

#options
my %options=();
getopts("abc",\%options);

#vary e-value cutoff (optional)
my $e_value=1;
if($options{c}){
    $e_value=$options{c};
}


#temporary directory
unless ($tmpdir=~/\/$/) {$tmpdir .="/";}
my $blast_dir_temp=$tmpdir;
if(!(-e $blast_dir_temp)){my $ret1=`mkdir $blast_dir_temp`;}

#copy subject file
print STDERR "copying subject file\n";
my $ret2=`cp $file_subject $blast_dir_temp/file_subject`;

#format subject file
print STDERR "formatting subject file\n";
my $ret3=`formatdb -i $blast_dir_temp/file_subject -p F`;

#blast query against database
print STDERR "blasting query file against subject file\n";
my $ret4=`blastall -p blastn -a 6 -d $blast_dir_temp/file_subject -i $file_query -o $blast_dir_temp/blastout -F F`;

#parse into blast_parsed format
print STDERR "parsing blast output\n";
my $ret5=`perl $perldir/blastoutparse.pl $blast_dir_temp/blastout > $blast_dir_temp/blastparsed`;
my $ret6=`perl $perldir/blastparselect.pl $blast_dir_temp/blastparsed -e > $blast_dir_temp/blastparsed_fullcov`;

#select output given the arguments
my $ret9;
unless($options{a} or $options{b}){$ret9=`cat $blast_dir_temp/blastparsed`;}
if($options{a}){$ret9=`cat $blast_dir_temp/blastparsed_fullcov`;}
if($options{b}){
    #sort according to miRDeep
    my $ret7=`cat $blast_dir_temp/blastparsed_fullcov | sort -n -k 6 > $blast_dir_temp/sorted`;
    my $ret8=`cat $blast_dir_temp/sorted | sort -k 4 > $blast_dir_temp/distribution`;
    $ret9=`cat $blast_dir_temp/distribution`;
}

#delete temporary directory
my $ret10=`rm -r $blast_dir_temp`;
my $ret11=`rm formatdb.log`;
#my $ret12=`rm -rf $tmpdir`;
#print
print $ret9;

exit;
