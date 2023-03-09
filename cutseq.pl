#!/usr/bin/perl -w

#           #
# cutseq.pl #
#           #


#ccccccccccccccccccccccccccccccccccccccccccccccc#
# Ran Friedman                                  #
# Computaional Chemistry and Biochemistry Group #
# Department of Chemistry and Biomedicinal      #
# Sciences, Linnaeus University, Kalmar, Sweden #
# email: ran.friedman@lnu.se                    #
#ccccccccccccccccccccccccccccccccccccccccccccccc#

use strict;
use Getopt::Long;
use Bio::SeqIO;

#
# Declaration and init
#

my $prog="cutseq.pl";
my $ver=0.1;
my $prog_install=".";
my $author_info="Ran Friedman\nComputational Chemistry and Biochemistry Group\nDepartment of Chemistry and Biomedicinal Sciences\nLinnaeus University\n";
my $flags=
"  -h/--help print help message and quit\n".
"  -s/--seqfile [] input sequence\n".
"  --format/--seqformat format of the input sequence. Default is to guess by \n".
"      extension or content and use fasta otherwise\n".
"  -b [i] where to begin \n".
"  -e [i] where to stop \n".
"  -o/--out [] output sequence\n".
"  --outformat [fasta] output format\n".
"";
my $usage="$prog\nUsage:\n  $prog (flags)\n$flags\n";
my ($help,$seqfile,$seqformat,$outfile,$outformat,$begin,$end,$bDEBUG);
my ($seq_obj,$seqio_obj,$seq_out,$seqstr);
		    
#
# main part of the program
#

if (@ARGV==0) {
  print $usage;
  exit (0);
}

my $opts= GetOptions( "help|h" => \$help,
		      "seqfile|s=s" => \$seqfile,
		      "seqformat|format=s" => \$seqformat,
		      "b=i" => \$begin,
		      "e=i" => \$end,
		      "out=s" => \$outfile,
		      "outformat=s" => \$outformat,
		      "d|debug" => \$bDEBUG
		    );

check_input();

if ($seqformat) {
  $seqio_obj = Bio::SeqIO->new(-file => $seqfile, 
			       -format => $seqformat);
}
else {
  $seqio_obj = Bio::SeqIO->new(-file => $seqfile);
}

$seq_obj = $seqio_obj->next_seq;
if ($begin && $end) { 
  $seqstr   = $seq_obj->subseq($begin,$end);
}
elsif ($begin) {
  $end=$seq_obj->length;
  $seqstr   = $seq_obj->subseq($begin,$end);
}
elsif ($end) {
  
   $seqstr   = $seq_obj->subseq(1,$end);
}
else {
  die "Cannot deal with the input: -b and -e values missing\n";
}

$seq_obj = Bio::PrimarySeq->new( -seq => $seqstr);
$seq_out = Bio::SeqIO->new(
			   -file   => ">$outfile",
			   -format => $outformat,
			  );
$seq_out->write_seq($seq_obj);

#
# subroutines
#

sub check_input {
  my $parameters="The program will use the following parameters:\n";
  
  if ($help) {
    print_help_message();
    exit(0);	
  }
  unless ($seqfile) {
    print STDERR "$prog error: no sequence file\nAborting.\n";
    exit(6);
  }
  $parameters.="\t--seqfile=$seqfile\n";
  if ($seqformat) {
    $parameters.="\t--seqformat=$seqformat\n";
  }
  if ($begin) {
    $parameters.="\t-b=$begin\n";
  }
  if ($end) {
     $parameters.="\t-b=$end\n";
  }
  unless ($outformat) {
    $outformat="fasta";
  }
  $parameters.="\t-outformat=$outformat\n";
  print $parameters,"\n";
}

sub print_help_message {
    my $help=
"The program $prog reads a protein or nucleic acid sequence and cuts it, \n".
"producing a sub-sequence that begins and ends at sub points given by -b and -e\n";
    print "$usage\n";
    print "$help\n";
    print "Programmed by:\n";
    print "$author_info\n";
}

0;
