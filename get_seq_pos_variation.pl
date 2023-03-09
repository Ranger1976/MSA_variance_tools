#!/usr/bin/perl -w

#                          #
# get_seq_pos_variation.pl #
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

my $prog="get_seq_pos_variation.pl";
my $ver=0.1;
my $prog_install=".";
my $author_info="Ran Friedman\nComputational Chemistry and Biochemistry Group\nDepartment of Chemistry and Biomedicinal Sciences\nLinnaeus University\n";
my $flags=
"  -h/--help print help message and quit\n".
"  --aln/--alnfile [] input (mRNA) sequence\n".
"  --format/--alnformat format of the input MSA. Default is to guess by \n".
"      extension \n".
"  --shift/shift-res [i] shift the first residue of the first sequence. This \n".
"      is useful when the alignment was carried out for a subseq of the first \n".
"      sequence, but mutation or residue in the input are given accrding to\n".
"      the canonical numbering. For example, if you carried out the MSA using\n".
"      only residues 38-147 of your protein you need to set shift-res to 38\n".
"  --res/--residues [xx yy zz ...] calculate variance only for these residues\n".
"  --resfile [] same as -res but take the residues from a file\n".
"    the file shall contain a list of whitespace separated residues.\n".
"    if --resfile is given --res will be ignored \n".
"  --mut/--mutfile [] for each mutation in the file, check if the variation\n".
"    is located in the MSA\n".
"    The format of the mutation file is:\n".
"    Number AA NN\n".
"    where AA is the original amino acid and NN is the new one\n".
"    For example: \n".
"    247 K N\n".
"    if --mut is given --res and --resfile will be ignored \n".  
"   -v Be verbose (more info printed)\n".
"";
my $usage="$prog\nUsage:\n  $prog (flags)\n$flags\n";
my ($help,$alnfile,$alnformat,$shiftres,$bDEBUG,$bMut,$bVerbose);
my @res_list;
my %mut_list;
#my ($mut_list_size,$mut_list_match)=(0,0);
my ($resfile,$mutfile);
my $list_of_residues;
my ($aln,$aln_io,$curres,$seq,$curseq,$res,$sum,$seqname,$slen);
my ($i,$j,$pos,$resnum,$nseq,$orig_res,$firstseq,$resi);
my ($cmutations_r,@cmutations);
my %count;
#my @pos; # pos[i] is the position of residue [i+1] in the MSA
		    
#
# main part of the program
#

if (@ARGV==0) {
  print $usage;
  exit (0);
}

my $opts= GetOptions( "help|h" => \$help,
		      "alnfile|s=s" => \$alnfile,
		      "alnformat|format=s" => \$alnformat,
		      "shift|shift-res=i" => \$shiftres,
		      "res|residues:i{,}" => \@res_list,
		      "resfile:s" => \$resfile,
		      "mut|mutfile=s" => \$mutfile,
		      "v|verbose" => \$bVerbose,
		      "d|debug" => \$bDEBUG
		    );

check_input();

if ($alnformat) {
  $aln_io = Bio::AlignIO->new(-file   => $alnfile,
			      -format => $alnformat);
}
else {
  # guess the format and hope for the best
  $aln_io = Bio::AlignIO->new(-file   => $alnfile);
}

$aln = $aln_io->next_aln();
$curseq = $aln->get_seq_by_pos(1);
$seqname = $curseq->id;
$firstseq=$curseq->seq;
$slen=$curseq->end;

#for ($i=0;$i<$slen;$i++) {
#  # generate the array of positions just once
#  push(@pos,$aln->column_from_residue_number($seqname, $i));
#}

# Work over each sequence and position
#foreach $seq ($aln->each_seq) {
if(%mut_list) {
  # have to go over the residues and mutations and do something with them
  foreach $resi (sort {$a<=>$b} keys %mut_list) {
    $cmutations_r=$mut_list{$resi};
    @cmutations=@$cmutations_r;
    # get the right residue
    $i=$resi-$shiftres+1;
    $pos=$aln->column_from_residue_number($seqname, $i);
    $orig_res=substr($firstseq,$pos-1,1);
    for ($j=0;$j<$#cmutations;$j+=2) { 
      if ($orig_res ne $cmutations[$j]) {
	print STDERR "$prog: Error! residue number $resi is ($orig_res and not $cmutations[$j] \n";
      }
      else { 
	# work over all the sequences
	# see if we have this mutation
	%count=();
	$nseq=0; # number of variations
	$bMut=0;
	foreach $seq ($aln->each_seq) {
	  $res = $seq->subseq($pos, $pos);
	  if ($res =~ /[ACDEFGHIKLMNPQRSTVWY]/) {
	    $count{$res}++;
	    $nseq++;
	  }
	  if ($res eq $cmutations[$j+1]) {
	    $bMut=1;
	    if ($bVerbose) {
	      printf "Sequence %s carried out the mutation %d %s -> %s\n",$seq->id,$resi,$orig_res,$cmutations[$j+1];
	    }
	  }
	}
	if ($bMut) {
	  printf "Mutation %d %s -> %s is present in (%.2f%%) of the sequences\n",$resi,$orig_res,$cmutations[$j+1], $count{$cmutations[$j+1]}/$nseq*100;
	}
	else {
	  printf "Mutation %d %s  -> %s is not present\n",$resi,$orig_res,$cmutations[$j+1];
	}
      }
    }
  }
}
elsif(@res_list) {
  # have to go over the residues 
  while ($resi=shift(@res_list)) {
    $i=$resi-$shiftres+1;
    ### This code should move to a sub ###
    $pos=$aln->column_from_residue_number($seqname, $i);
    $orig_res=substr($firstseq,$pos-1,1);
    $resnum=$i+$shiftres-1;
    %count=();
    $nseq=0; # number of variations
    foreach $seq ($aln->each_seq) {
      $res = $seq->subseq($pos, $pos);
      # Now we deal with stats only for real residues
      if ($res =~ /[ACDEFGHIKLMNPQRSTVWY]/) {
	$count{$res}++;
	$nseq++;
      }
    }
    printf "%s%d variations ",$orig_res,$resnum;
    foreach $res (keys %count) {
      printf "%s %d (%d%%) ",$res, $count{$res}, $count{$res}/$nseq*100;
    }
    print "\n";
  }
  ### Up till here ###
}
else {
  for ($i=1;$i<=$slen;$i++) {
    ### This code should move to a sub ###
    $pos=$aln->column_from_residue_number($seqname, $i);
    $orig_res=substr($firstseq,$pos-1,1);
    $resnum=$i+$shiftres-1;
    %count=();
    $nseq=0; # number of variations
    foreach $seq ($aln->each_seq) {
      $res = $seq->subseq($pos, $pos);
      # Now we deal with stats only for real residues
      if ($res =~ /[ACDEFGHIKLMNPQRSTVWY]/) {
	$count{$res}++;
	$nseq++;
      }
    }
    printf "%s%d variations ",$orig_res,$resnum;
    foreach $res (keys %count) {
      printf "%s %d (%d%%) ",$res, $count{$res}, $count{$res}/$nseq*100;
    }
    print "\n";
    ### Up till here ###
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
  unless ($shiftres) {
    $shiftres=1;
  }
  $parameters.="\t--shiftres=$shiftres (1 means no shift)\n";
  #unless ($n_subst) {
  #  $n_subst=1;
  #}
  #$parameters.="\t--number-of-substitutions=$n_subst\n";
  if ($mutfile) {
    undef @res_list;
    $parameters.="\t--mutfile =$mutfile\n";
    %mut_list=process_mutations_file($mutfile);
  }
  elsif ($resfile) {
    undef @res_list;
    $parameters.="\t--resfile =$resfile\n";
    $parameters.="\t parameters in --res/--residues will be ignored (if given)\n";
    @res_list=process_residue_file($resfile);
    if ($bDEBUG) {
      $list_of_residues=ar2str(\@res_list);
      print "Residues to modify $list_of_residues\n";
    }
  }
  else {
    if (@res_list) {
      $list_of_residues=ar2str(\@res_list);
      $parameters.="\t--residues =$list_of_residues\n";
    }
  }
  print $parameters,"\n";
}

sub print_help_message {
    my $help=
"The program $prog get a multiple sequence alignment as an input, and optionally\n".
"also a list of potential mutations in a file. \n".
"If a list of mutations is given, the program will search the MSA for each \n".
"mutation, to check if this mutation is manifested at the MSA and how common it\n".
"is. Otherwise, the program lists all possible variations for the first sequence\n".
"in the MSA or for a list of residues given in a separate file or command line\n".
"Examples:\n".
"";
#"See the manual for a more details\n";
    print "$usage\n";
    print "$help\n";
    print "Programmed by:\n";
    print "$author_info\n";
}

sub process_residue_file {
  my @ret;
  my $resfile=shift;
  my @ar;
  my $res;
  open (RES,$resfile)||die "$prog: Error! cannot open file $resfile :$!\n";
  while ($_ = <RES>) {
    next if (/^[#;]/);
    @ar=split;
    while ($res=shift(@ar)) {
      if ($res=~/^\d+$/) {
	push(@ret,$res);
      }
    }
  }
  close (RES);
  return @ret;
}

sub process_mutations_file {
  # The output is a hash:
  # residue -> \([AA] [BB] )
  # if a residue is involved in more than one mutation
  # residue -> \([AA] [BB] [AA] [CC] ...)
  my %ret;
  my $resfile=shift;
  my $res;
  open (RES,$resfile)||die "$prog: Error! cannot open file $resfile :$!\n";
  while ($_ = <RES>) {
    next if (/^[#;]/);
    if (/(\d+)\s+([A-Z])\s+([A-Z])/) {
      my @ar=($2,$3);
      if (exists $ret{$1}) {
	my $ar2_r=$ret{$1};
	my @ar2=@$ar2_r;
	unshift (@ar,@ar2);
      }
      $ret{$1}=\@ar;
    }
  }
  close (RES);
  return %ret;
}

sub ar2str {
  #input: reference to an array, separator (optional)
  #returns: a string with all values
  # example: 
  # input: ar2str(\[1,2,3], " , ") 
  # output: 1 , 2 , 3
  my $ar_r=shift;
  my $sep = @_ ? shift : " ";
  my @ar=@$ar_r;
  my $i;
  my $ret="";
  for ($i=0;$i<$#ar;$i++) {
    $ret.=$ar[$i].$sep;
  }
  $ret.=$ar[$#ar];
  return $ret;
}

0;
