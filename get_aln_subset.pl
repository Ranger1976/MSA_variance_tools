#!/usr/bin/perl -w

#                          #
# get_aln_subset.pl #
#                          #

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
use Bio::AlignIO;

#
# Declaration and init
#

my $prog="get_aln_subset.pl";
my $ver=0.1;
my $prog_install=".";
my $author_info="Ran Friedman\nComputational Chemistry and Biochemistry Group\nDepartment of Chemistry and Biomedicinal Sciences\nLinnaeus University\n";
my $flags=
"  -h/--help print help message and quit\n".
"  --aln/--alnfile [s] input (protein) sequence alignment file\n".
"  --format format of the input MSA. Default is to guess by \n".
"      extension \n".
"  --first/--firstres [i] first residue to consider \n".
"  --last/--lastres [i] last residue to consider \n".
"  --lnt/--minlength [i] minimum length aligned (bp) \n".
"  --pct/--minlengthpct [f] minimum length aligned (pct) \n".
"    if both --lnt and --lpct are given, --lnt is ignored \n".
"  --o/--out [s] output file\n".
"    if not such file is provided the sequences will be written to the \n".
"    standard output\n".
"   -v Be verbose (more info printed)\n".
"\n".
"Note that the reference point is always the first sequence.\n".
"";
my $usage="$prog\nUsage:\n  $prog (flags)\n$flags\n";
my ($help,$alnfile,$alnformat,$firstres,$lastres,$minlnt,$bDEBUG,$bVerbose);
my ($lntpct,$outfile,$outseq);
my ($fh,$bFILEOPEN);
my ($aln,$aln_io,$seq,$curseq,$seqname,$slen,$res,$firstseq);
my $seqaln_lnt;
my $i;
my @col;
my $pos;
my @approved_seq;
		    
#
# main part of the program
#

if (@ARGV==0) {
  print $usage;
  exit (0);
}

my $opts= GetOptions( "help|h" => \$help,
		      "alnfile|aln=s" => \$alnfile,
		      "format=s" => \$alnformat,
		      "firstres|first=i" => \$firstres,
		      "lastres|last=i" => \$lastres, 
		      "minlength|lnt=i" => \$minlnt,
		      "minlengthpct|pct=f"  => \$lntpct,
		      "out|o=s" => \$outfile,
		      "v|verbose" => \$bVerbose,
		      "d|debug" => \$bDEBUG
		    );

check_input();

# Generate the aln object

if ($alnformat) {
  $aln_io = Bio::AlignIO->new(-file   => $alnfile,
			      -format => $alnformat);
}
else {
  # guess the format and hope for the best
  $aln_io = Bio::AlignIO->new(-file   => $alnfile);
}

# Get the first sequence

$aln = $aln_io->next_aln();
$curseq = $aln->get_seq_by_pos(1);
$seqname = $curseq->id;
$firstseq=$curseq->seq;
$slen=$curseq->end;

print "Reference sequence $seqname length $slen\n" if $bVerbose;

if ($minlnt && $minlnt > $slen) {
    print STDERR "$0: Error in input: the minimum requested length $minlnt is larger\n". 
	"than the length of the first sequence with is $slen\n";	
    print STDERR "Aborting $0";
    exit(6);
}

if  ($lastres && $lastres > $slen) {
    print STDERR "$0: Error in input: the last position $lastres is not found in the first \n".
	"sequence";
    print STDERR "Aborting $0";
    exit(6);
}

if ($lntpct) {
    if ($firstres && $lastres) {
	$minlnt=($lastres-$firstres+1)*$lntpct/100;
    }
    elsif ($firstres) {
	$minlnt=($slen-$firstres+1)*$lntpct/100;
    }
    elsif ($lastres) {
	$minlnt=$lastres*$lntpct/100;
    }
    else {
	$minlnt=$slen*$lntpct/100;
    }
}
$minlnt = defined($minlnt) ? $minlnt : 1;
if ($bDEBUG) {
    print "minlength is $minlnt \n";
}

# After the first seq is read, we will generate an array of columns that includes this sequence or a subsequence between $firstres and $lastres
$firstres = $firstres ? $firstres : 1;
$lastres = $lastres ? $lastres: $slen;
for ($i=$firstres ; $i<=$lastres; $i++) {
    push(@col,$aln->column_from_residue_number($seqname, $i));
}

# Go over each sequence and check alignment length to the first
foreach $curseq ($aln->each_seq) {
    $seqaln_lnt=0;
    foreach $pos (@col) {
	$res = $curseq->subseq($pos, $pos);
	if ($res =~ /[ACDEFGHIKLMNPQRSTVWY]/) {
	    $seqaln_lnt++;
	}
    }
    if ($seqaln_lnt >= $minlnt) {
	#push (@approved_seq,$curseq->id);
	print $curseq->display_id, " is aligned\n" if $bVerbose;
	$outseq->write_seq($curseq);
    }
}

#
# subroutines
#

sub check_input {
  my $parameters="The program will use the following parameters:\n";
  
  if ($help) {
    print_help_message();
    exit(0);	
  }
  unless ($alnfile) {
    print STDERR "$prog error: no sequence file\nAborting.\n";
    exit(6);
  }
  $parameters.="\t--alnfile=$alnfile\n";
  if ($alnformat) {
    $parameters.="\t--alnformat=$alnformat\n";
  }
  if ($firstres){
      if ($firstres<1) {
	  print STDERR "First residue cannot be $firstres\n";
	  print STDERR "Ignoring parameter --firstres\n";
	  $firstres=0;
      }
      else {
	  $parameters.="\t--firstres=$firstres\n";
      }
  }
  if ($lastres) {
      if ($lastres<1) {
	  print STDERR "Last residue cannot be $lastres\n";
	  print STDERR "Ignoring parameter --lastres\n";
	  $lastres=0;
      }
      elsif ($lastres<$firstres) {
	  print STDERR "The last residue cannot be smaller than the first \n";
	  print STDERR "Ignoring parameter --lastres\n";
	  $lastres=0;
      }
      else {
	  $parameters.="\t--lastres=$lastres\n";
      }
  }
  if ($minlnt) {
      if ($minlnt<1) {
	  print STDERR "Sequence length cannot be $minlnt\n";
	  print STDERR "Ignoring parameter --minlength\n";
	  $minlnt=0;
      }
      elsif ($firstres && $lastres && $minlnt > $lastres-$firstres+1) {
	  print STDERR "The minimum length of the alignment cannot be larger than the size of the input sequence\n";
	  print STDERR "Ignoring parameter --minlength\n";
	  $minlnt=0;
      }
      else {
	  $parameters.="\t--minlength=$minlnt\n";
      }
  }
  if ($lntpct) {
      if ($lntpct<0 || $lntpct>100) {
	  print STDERR "The percentage of the aligned sequence must be between 0 and 100\n";
	  $lntpct=0;
      }
      else {
	  $parameters.="\t--minlengthpct=$lntpct\n";
      }
  }
  unless ($firstres || $lastres || $lntpct || $minlnt) {
      print STDERR "Error in input - nothing to do\n";
      print STDERR "You must supply at list --firstres, --lastres, --minlength or --minlengthpct\n";
      print STDERR "Program $0 aborts\n";
      exit (6);
  }
  if ($minlnt && $lntpct && $bVerbose) {
      print "Both --minlength and --minlengthpct were given \n";
      print "\t--minlengthpct will be used";
      print "\t--minlength will be ignored";
  }
  if ($outfile) {
      #$bFILEOPEN=open ($fh, ">",$outfile); 
      #unless ($bFILEOPEN) {
	#  print STDERR "Error in input - cannot open file $outfile\n";
	#  print STDERR "$!\n";
	#  print STDERR "Program $0 aborts\n";
	#  exit (6);
      #}
      #close ($fh);
      $outseq=Bio::SeqIO->new( -file   => ">$outfile",
			       -format => "fasta",
	  );
  }
  else {
      $outseq = Bio::SeqIO->new( -fh     => \*STDOUT,
				 -format => "fasta",
	  );
  }
  print $parameters,"\n";
}

sub print_help_message {
    $help=
"The program $prog get a multiple sequence alignment as an input, and shall \n".
"output a set sequences the fullfill a certain coverage criterion.\n".
"This could be:\n".
"1. Sequences that have their alignment from a certain residue AA up to a, \n".
"   certain residue BB, or between AA and BB.\n".
"2. Have an alignment of certain length LL or more.\n".
"3. Have an alignment that covers at least PP % of the length.\n".
"4. A combination of (1) and (2) or  (1) and (3), e.g., an alignment that \n".
"   covers at least PP % of the sequence between positions AA and BB\n". 
"   \n".
"Examples:\n".
"$0 --first 100 --pct 50 --aln my.aln -o out.fa \n".
"   all sequences that cover 50% of the area starting at position 100 of the\n"."   first sequence\n".
"$0 --first 50 --last 550 --length 400 --aln my.aln -o out.fa \n".
"   all sequences that have an alignment of at least 400 between 50 and 550\n".
"   of the first sequence\n".
"$0 --first 50 --last 550 --aln my.aln -o out.fa \n".
"   all sequences that are aligned between positions 50 and 550 of the first\n".
"   sequence\n".
"\n".
"The alignment is assumed to contain protein sequences, but this is not \n".
"verified.\n".
"The output is a fasta format file\n";
    print "$usage\n";
    print "$help\n";
    print "Programmed by:\n";
    print "$author_info\n";
}

0;
