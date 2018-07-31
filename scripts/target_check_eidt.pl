#!/usr/bin/perl -w
#Filename:
#Author: Tian Dongmei
#Email: tiandm@big.ac.cn
#Date: 2009-05-06
#Modified:
#Description: É¾³ýmatched reads
my $version = 1.00;

use strict;
use Getopt::Long;

my ($infile) = @ARGV;
my ( %opts, $outdir );
open( OPT, "<$infile" ) or die "cannot open config.properties\n";

while ( my $line = <OPT> ) {
    if ( $line =~ /RESULT_DIRECTORY\s*=\s*(\S+)/ ) { $opts{'out'} = $1; }
}
close OPT;
$opts{'mi'} = $opts{'out'} . "\/target\/miranda_targets";
$opts{'hy'} = $opts{'out'} . "\/target\/RNAhybrid_targets";

GetOptions( \%opts, "p=s" );

if ( !( defined $opts{mi} and defined $opts{hy} ) ) {    #necessary arguments
    &usage;
}

my $filein  = $opts{'mi'};
my $filein2 = $opts{'hy'};
my $fileout = $opts{'out'} . "\/target\/target_predict";

my %uniq;

my %hash;
open IN, "<$filein";                                     #input file
while ( my $aline = <IN> ) {
    chomp $aline;
    my @tmp = split /\t/, $aline;
    $tmp[0] =~ s/>>//;
    my @tar = split /\|/, $tmp[1];
    $hash{ $tmp[0] }->{ $tar[0] } = $tmp[4];
}
close IN;

my %target;
open IN, "<$filein2";
while ( my $aline = <IN> ) {
    chomp $aline;
    if ( $aline =~ /^>/ ) {
        my @tmp = split /\t/, $aline;
        $tmp[0] =~ s/>//;
        my @tar = split /\|/, $tmp[1];
        if ( defined $hash{ $tmp[0] }->{ $tar[0] } ) {

            #	push @{$target{$tmp[0]}},$tar[0]};
            $uniq{ $tmp[0] }->{ $tar[0] } = $hash{ $tmp[0] }->{ $tar[0] };
            delete( $hash{ $tmp[0] }->{ $tar[0] } );
        }
    }
}
close IN;

open OUT, ">$fileout";

#foreach my $key (keys %target) {
#	print OUT "$key\t";
#	print OUT scalar @{$target{$key}};
#	print OUT "\t@{$target{$key}}\n";
#}
foreach my $key ( keys %uniq ) {
    print OUT "$key\t";
    my @value =
      sort { $uniq{$key}{$b} <=> $uniq{$key}{$a} } keys %{ $uniq{$key} };
    print OUT scalar @value;
    print OUT "\t@value\n";
}
close OUT;

sub usage {
    print <<"USAGE";
Version $version
Usage:
$0 -mi -hy -o
options:
-mi input file,miRanda target predict result
-hy input file,RNAhybrid target predict result.
-o output file
-h help
USAGE
    exit(1);
}

