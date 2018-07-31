#!/usr/bin/perl 
#Author: Gong Jing 2011-12-6
# used to get miRNA precursor.txt AND miRNA.txt for miREXpress;
use Getopt::Long;

my %opts;
GetOptions( \%opts, "p=s", "m=s", "d=s", "g=s", "species", "h" );
if (
    !(
            defined $opts{p}		#p: hairpin.fa;
        and defined $opts{m}		#m: mature.fa;
        and defined $opts{d}  		#d: miRNA.dat;
		and defined $opts{g}		#g: hsa.gff; or other species gff; 
	#  and defined $opts{o}		#o:	hairpin.txt for miRExpress;
	#   and defined $opts{o2}		#o2:mature.txt for miRExpress;
	#	and defined $opts{o3}		#o3: miRNA.info
	)
    || $opts{h}
  )
{
    &usage;
}
open( OUT,  ">hairpin.txt" );	#plese mv the three files to the data file for following analysis;
open( OUT2, ">mature.txt" );
open( OUT3, ">miRNA.info");
if ( $opts{species} eq "" ) {
    $species = "hsa";
}
else {
    $species = $opts{species};
}

#get outfile1;
&parse_file_fasta( $opts{p}, "pre" );
foreach $mir ( keys %pre ) {
    if ( $mir =~ /$species/ ) {
        print OUT "$mir\t$pre{$mir}\n";
    }
}

#get outfile2;
&parse_file_fasta( $opts{m}, "mat" );

open( IN3, "<$opts{d}" );
$/ = "//";
while ( $record = <IN3> ) {
	if ( $record =~ /product\=.*product\=/s ) {

        $record =~
/ID\s+(\S+).*miRNA\s+(\d+)\.\.(\d+).*product="(\S+)".*miRNA\s+(\d+)\.\.(\d+).*product="(\S+)"/s;
       
my 	    $pre    = $1;
my      $mat1_s = $2;
my 		$mat1_e = $3;
my      $mat1   = $4;
my      $mat2_s = $5;
my		$mat2_e = $6;
my      $mat2   = $7;
        $hash{$mat1} = $pre;
        $hash{$mat2} = $pre;
        $pos{$mat1}{s}  = $mat1_s;
        $pos{$mat2}{s}  = $mat2_s;
		$pos{$mat1}{e}  = $mat1_e;
		$pos{$mat2}{e}  = $mat2_e;
	}
    else {
        $record =~ /ID\s+(\S+).*miRNA\s+(\d+)\.\.(\d+).*product="(\S+)"/s;
  my      $pre         = $1;
  my      $mat1_s      = $2;
  my   	  $mat1_e      = $3;
  my 	  $mat1        = $4;
          $hash{$mat1} = $pre;
          $pos{$mat1}{s}  = $mat1_s;
    	  $pos{$mat1}{e}  = $mat1_e
	  }
}
close IN3;
$/="\n";

open (IN4,"<$opts{g}") or die "cannot open $opts{g}\n";
my %info;
while(<IN4>){
	chomp;
	if(/#/){next;}
	my ($chr,$s,$e,$strand,$ID)=(split/\t/)[0,3,4,6,8];
	$chr="chr".$chr;
	$ID=~/ACC="(\S+)".*ID="(\S+)"/;
	my $acc=$1;
	my $id=$2;

	$descri="$chr\t$strand\t$s\t$e";
	$info{$id}=$descri;

}
close IN4;

print OUT3 "#mature_id\tpre_id\tchr\tstrand\tpre_start\tpre_end\tmat_start\tmature_end\tmat_seq\tpre_seq\n";
foreach $mat_mir ( keys %hash ) {
    if($mat_mir=~/$species/){
	print OUT2 "$mat_mir\t$mat{$mat_mir}\t$hash{$mat_mir}\t$pos{$mat_mir}{s}\n";
	my $pre_id=$hash{$mat_mir};
	print OUT3 "$mat_mir\t$pre_id\t$info{$pre_id}\t$pos{$mat_mir}{s}\t$pos{$mat_mir}{e}\t$mat{$mat_mir}\t$pre{$pre_id}\n";
}
}
close OUT3;


sub parse_file_fasta {
    my ( $file, $hash ) = @_;
    my ( $id, $desc, $sequence ) = ();

    open( FASTA, "<$file" ) or die "can not open $file\n";
    while (<FASTA>) {
        chomp;
        if (/^>(\S+)(.*)/) {
            $id       = $1;
            $desc     = $2;
            $sequence = "";
            while (<FASTA>) {
                chomp;
                if (/^>(\S+)(.*)/) {
                    $$hash{$id} = $sequence;
                    $id         = $1;
                    $desc       = $2;
                    $sequence   = "";
                    next;
                }
                $sequence .= $_;
            }
        }
    }
    $$hash{$id} = $sequence;
    close FASTA;
}

sub usage {
    print <<"USAGE";
	Usage:
	$0  -p -m -species -d -g 
	options:
	-p      precurors file
	-m      mature miRNA file
	-d		miRNA.dat
	-g		hsa.gff
	-species species alias name :such as hsa; default is hsa;
	 output files: precurors.txt  mature.txt miRNA.info 
	-h      help
USAGE
    exit(1);
}

