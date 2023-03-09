#!/usr/bin/perl -w

#          #
# sland.pl #
#          #

# SLand : Substitution LANDscape

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

my $prog="sland";
my $ver=0.1;
my $prog_install=".";
my $author_info="Ran Friedman\nComputational Chemistry and Biochemistry Group\nDepartment of Chemistry and Biomedicinal Sciences\nLinnaeus University\n";
my $flags=
"  -h/--help print help message and quit\n".
"  -s/--seqfile [] input (mRNA) sequence\n".
"  --format/--seqformat format of the input sequence. Default is to guess by \n".
"      extension or content and use fasta otherwise\n".
"  --first-aa [1] how should the first amino acid be numbered\n".
#"  -n/--number-of-substitutions [1] \n".
"  --res/--residues [xx yy zz ...] substitute only these residues\n".
"  --resfile [] same as -res but take the residues from a file\n".
"    the file shall contain a list of whitespace separated residues.\n".
"    if --resfile is given --res will be ignored \n".
"  --mut/--mutfile [] for each mutation in the file, check if it is\n".
"    a single nucleotid mutation\n".
"    The format of the mutation file is:\n".
"    Number AA NN\n".
"    where AA is the original amino acid and NN is the new one\n".
"    For example: \n".
"    247 K N\n".
"    if --mut is given --res and --resfile will be ignored \n".  
"  --exclude-stop ignore nonsense mutations\n".
"  --prep-mut-file prepare a file, based on the given mutations, Use of this \n".
"    flag makes --mut obsolete\n".
"";
my $usage="$prog\nUsage:\n  $prog (flags)\n$flags\n";
my ($help,$seqfile,$seqformat,$firstaa,$n_subst,$exclude_stop,$bDEBUG);
my ($prep_file);
my @res_list;
my %mut_list;
my ($mut_list_size,$mut_list_match)=(0,0);
my ($resfile,$mutfile);
my $list_of_residues;
my ($seq_obj,$seqio_obj,$prot_obj);
my ($orig_seq,$ppos,$plen,$codon,$res);
my $i;
my $cmutations_r;
my @cmutations;
my ($aa,@pos_subsitutions);


my %translation = (
		   'UCA' => 'S',    # Serine
		   'UCC' => 'S',    # Serine
		   'UCG' => 'S',    # Serine
		   'UCU' => 'S',    # Serine
		   'UUC' => 'F',    # Phenylalanine
		   'UUU' => 'F',    # Phenylalanine
		   'UUA' => 'L',    # Leucine
		   'UUG' => 'L',    # Leucine
		   'UAC' => 'Y',    # Tyrosine
		   'UAU' => 'Y',    # Tyrosine
		   'UAA' => '_',    # Stop
		   'UAG' => '_',    # Stop
		   'UGC' => 'C',    # Cysteine
		   'UGU' => 'C',    # Cysteine
		   'UGA' => '_',    # Stop
		   'UGG' => 'W',    # Tryptophan
		   'CUA' => 'L',    # Leucine
		   'CUC' => 'L',    # Leucine
		   'CUG' => 'L',    # Leucine
		   'CUU' => 'L',    # Leucine
		   'CCA' => 'P',    # Proline
		   'CCC' => 'P',    # Proline
		   'CCG' => 'P',    # Proline
		   'CCU' => 'P',    # Proline
		   'CAC' => 'H',    # Histidine
		   'CAU' => 'H',    # Histidine
		   'CAA' => 'Q',    # Glutamine
		   'CAG' => 'Q',    # Glutamine
		   'CGA' => 'R',    # Arginine
		   'CGC' => 'R',    # Arginine
		   'CGG' => 'R',    # Arginine
		   'CGU' => 'R',    # Arginine
		   'AUA' => 'I',    # Isoleucine
		   'AUC' => 'I',    # Isoleucine
		   'AUU' => 'I',    # Isoleucine
		   'AUG' => 'M',    # Methionine
		   'ACA' => 'T',    # Threonine
		   'ACC' => 'T',    # Threonine
		   'ACG' => 'T',    # Threonine
		   'ACU' => 'T',    # Threonine
		   'AAC' => 'N',    # Asparagine
		   'AAU' => 'N',    # Asparagine
		   'AAA' => 'K',    # Lysine
		   'AAG' => 'K',    # Lysine
		   'AGC' => 'S',    # Serine
		   'AGU' => 'S',    # Serine
		   'AGA' => 'R',    # Arginine
		   'AGG' => 'R',    # Arginine
		   'GUA' => 'V',    # Valine
		   'GUC' => 'V',    # Valine
		   'GUG' => 'V',    # Valine
		   'GUU' => 'V',    # Valine
		   'GCA' => 'A',    # Alanine
		   'GCC' => 'A',    # Alanine
		   'GCG' => 'A',    # Alanine
		   'GCU' => 'A',    # Alanine
		   'GAC' => 'D',    # Aspartic Acid
		   'GAU' => 'D',    # Aspartic Acid
		   'GAA' => 'E',    # Glutamic Acid
		   'GAG' => 'E',    # Glutamic Acid
		   'GGA' => 'G',    # Glycine
		   'GGC' => 'G',    # Glycine
		   'GGG' => 'G',    # Glycine
		   'GGU' => 'G'     # Glycine
		  );

my %substitution = (
		    'AAA' => 'E I N Q R T _ ',
		    'AAC' => 'D H I K S T Y ',
		    'AAG' => 'E M N Q R T _ ',
		    'AAU' => 'D H I K S T Y ',
		    'ACA' => 'A I K P R S ',
		    'ACC' => 'A I N P S ',
		    'ACG' => 'A K M P R S ',
		    'ACU' => 'A I N P S ',
		    'AGA' => 'G I K S T _ ',
		    'AGC' => 'C G I N R T ',
		    'AGG' => 'G K M S T W ',
		    'AGU' => 'C G I N R T ',
		    'AUA' => 'K L M R T V ',
		    'AUC' => 'F L M N S T V ',
		    'AUG' => 'I K L R T V ',
		    'AUU' => 'F L M N S T V ',
		    'CAA' => 'E H K L P R _ ',
		    'CAC' => 'D L N P Q R Y ',
		    'CAG' => 'E H K L P R _ ',
		    'CAU' => 'D L N P Q R Y ',
		    'CCA' => 'A L Q R S T ',
		    'CCC' => 'A H L R S T ',
		    'CCG' => 'A L Q R S T ',
		    'CCU' => 'A H L R S T ',
		    'CGA' => 'G L P Q _ ',
		    'CGC' => 'C G H L P S ',
		    'CGG' => 'G L P Q W ',
		    'CGU' => 'C G H L P S ',
		    'CUA' => 'I P Q R V ',
		    'CUC' => 'F H I P R V ',
		    'CUG' => 'M P Q R V ',
		    'CUU' => 'F H I P R V ',
		    'GAA' => 'A D G K Q V _ ',
		    'GAC' => 'A E G H N V Y ',
		    'GAG' => 'A D G K Q V _ ',
		    'GAU' => 'A E G H N V Y ',
		    'GCA' => 'E G P S T V ',
		    'GCC' => 'D G P S T V ',
		    'GCG' => 'E G P S T V ',
		    'GCU' => 'D G P S T V ',
		    'GGA' => 'A E R V _ ',
		    'GGC' => 'A C D R S V ',
		    'GGG' => 'A E R V W ',
		    'GGU' => 'A C D R S V ',
		    'GUA' => 'A E G I L ',
		    'GUC' => 'A D F G I L ',
		    'GUG' => 'A E G L M ',
		    'GUU' => 'A D F G I L ',
		    'UAA' => 'E K L Q S Y ',
		    'UAC' => 'C D F H N S _ ',
		    'UAG' => 'E K L Q S W Y ',
		    'UAU' => 'C D F H N S _ ',
		    'UCA' => 'A L P T _ ',
		    'UCC' => 'A C F P T Y ',
		    'UCG' => 'A L P T W _ ',
		    'UCU' => 'A C F P T Y ',
		    'UGA' => 'C G L R S W ',
		    'UGC' => 'F G R S W Y _ ',
		    'UGG' => 'C G L R S _ ',
		    'UGU' => 'F G R S W Y _ ',
		    'UUA' => 'F I S V _ ',
		    'UUC' => 'C I L S V Y ',
		    'UUG' => 'F M S V W _ ',
		    'UUU' => 'C I L S V Y '
		   );

my %substitution_no_stop = (
			    'AAA' => 'E I N Q R T ',
			    'AAC' => 'D H I K S T Y ',
			    'AAG' => 'E M N Q R T ',
			    'AAU' => 'D H I K S T Y ',
			    'ACA' => 'A I K P R S ',
			    'ACC' => 'A I N P S ',
			    'ACG' => 'A K M P R S ',
			    'ACU' => 'A I N P S ',
			    'AGA' => 'G I K S T ',
			    'AGC' => 'C G I N R T ',
			    'AGG' => 'G K M S T W ',
			    'AGU' => 'C G I N R T ',
			    'AUA' => 'K L M R T V ',
			    'AUC' => 'F L M N S T V ',
			    'AUG' => 'I K L R T V ',
			    'AUU' => 'F L M N S T V ',
			    'CAA' => 'E H K L P R ',
			    'CAC' => 'D L N P Q R Y ',
			    'CAG' => 'E H K L P R ',
			    'CAU' => 'D L N P Q R Y ',
			    'CCA' => 'A L Q R S T ',
			    'CCC' => 'A H L R S T ',
			    'CCG' => 'A L Q R S T ',
			    'CCU' => 'A H L R S T ',
			    'CGA' => 'G L P Q ',
			    'CGC' => 'C G H L P S ',
			    'CGG' => 'G L P Q W ',
			    'CGU' => 'C G H L P S ',
			    'CUA' => 'I P Q R V ',
			    'CUC' => 'F H I P R V ',
			    'CUG' => 'M P Q R V ',
			    'CUU' => 'F H I P R V ',
			    'GAA' => 'A D G K Q V ',
			    'GAC' => 'A E G H N V Y ',
			    'GAG' => 'A D G K Q V ',
			    'GAU' => 'A E G H N V Y ',
			    'GCA' => 'E G P S T V ',
			    'GCC' => 'D G P S T V ',
			    'GCG' => 'E G P S T V ',
			    'GCU' => 'D G P S T V ',
			    'GGA' => 'A E R V ',
			    'GGC' => 'A C D R S V ',
			    'GGG' => 'A E R V W ',
			    'GGU' => 'A C D R S V ',
			    'GUA' => 'A E G I L ',
			    'GUC' => 'A D F G I L ',
			    'GUG' => 'A E G L M ',
			    'GUU' => 'A D F G I L ',
			    'UAA' => 'E K L Q S Y ',
			    'UAC' => 'C D F H N S ',
			    'UAG' => 'E K L Q S W Y ',
			    'UAU' => 'C D F H N S ',
			    'UCA' => 'A L P T ',
			    'UCC' => 'A C F P T Y ',
			    'UCG' => 'A L P T W ',
			    'UCU' => 'A C F P T Y ',
			    'UGA' => 'C G L R S W ',
			    'UGC' => 'F G R S W Y ',
			    'UGG' => 'C G L R S ',
			    'UGU' => 'F G R S W Y ',
			    'UUA' => 'F I S V ',
			    'UUC' => 'C I L S V Y ',
			    'UUG' => 'F M S V W ',
			    'UUU' => 'C I L S V Y '
			   );
		    
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
		      "first-aa=i" => \$firstaa,
		      #"number-of-substitutions|n=i"  => \$n_subst,
		      "res|residues:i{,}" => \@res_list,
		      "resfile:s" => \$resfile,
		      "mut|mutfile=s" => \$mutfile,
		      "exclude-stop" => \$exclude_stop,
		      "prep-mut-file=s" => \$prep_file,
		      "d|debug" => \$bDEBUG
		    );

check_input();

if ($prep_file) {
  open (MUT,">$prep_file") || die "Cannot open file $prep_file for writing: $!\n";
}

if ($seqformat) {
  $seqio_obj = Bio::SeqIO->new(-file => $seqfile, 
				   -format => $seqformat,
				   -alphabet => 'rna');
}
else {
  $seqio_obj = Bio::SeqIO->new(-file => $seqfile,
				  -alphabet => 'rna');
}

while ($seq_obj = $seqio_obj->next_seq){
  $orig_seq=$seq_obj->seq;
  $prot_obj = $seq_obj->translate;
  $plen=$prot_obj->length();
  if($bDEBUG) {
    print "Protein length (incl. stop if in sequence) $plen\n";
  }
  #print $prot_obj->seq."\n";
  if(%mut_list) {
    foreach $res (sort {$a<=>$b} keys %mut_list) {
      if ($res > $plen) {
	print STDERR "$prog: Error! residue number $res does not exist\n";
      }
      else {
	$codon=substr($orig_seq,($res-$firstaa)*3,3);
	$codon =~ tr/T/U/;
	$cmutations_r=$mut_list{$res};
	@cmutations=@$cmutations_r;
	for ($i=0;$i<$#cmutations;$i+=2) {
	  $mut_list_size++; # we count every occurrence
	  if ($translation{$codon} ne $cmutations[$i]) {
	    print STDERR "$prog: Error! residue number $res is $translation{$codon} and not $cmutations[$i] \n";
	  }
	  else {
	    # go through the list of possible mutations
	    # and check if the right residue exists there
	    if ($substitution{$codon} =~ /$cmutations[$i+1]/) {
	      print "Mutation ",$res," ",$cmutations[$i]," -> " ,$cmutations[$i+1]," is a single nucleotide substitution\n";
	      $mut_list_match++;
	    }
	    else {
	      print "Mutation ",$res," ",$cmutations[$i]," -> " ,$cmutations[$i+1]," is NOT a single nucleotide substitution\n";
	    }
	  }
	}
      }
    }
    print "\n $mut_list_match out of $mut_list_size mutations match\n";
  }
  elsif(@res_list) {
    # Deal only with the list of residues
    while ($res=shift(@res_list)) {
      last unless ($res < $plen);
      $codon=substr($orig_seq,($res-$firstaa)*3,3);
      $codon =~ tr/T/U/;
      if ($exclude_stop) {
	print $res," ",$translation{$codon}," -> ",$substitution_no_stop{$codon},"\n";
      }
      else {
	print $res," ",$translation{$codon}," -> ",$substitution{$codon},"\n";
      }
    }
  }
  else {
    for ($i=0;$i<$plen;$i++) {
      $ppos=$firstaa+$i;
      # go over the RNA seq, extract codons and check mutations
      $codon=substr($orig_seq,$i*3,3);
      $codon =~ tr/T/U/;
      if ($exclude_stop) {
	print $ppos," ",$translation{$codon}," -> ",$substitution_no_stop{$codon},"\n";
      }
      else {
	print $ppos," ",$translation{$codon}," -> ",$substitution{$codon},"\n";
      }      
      if ($prep_file) {
	$_= $exclude_stop ? $substitution_no_stop{$codon} : $substitution{$codon};
	@pos_subsitutions = split;
	foreach $aa (@pos_subsitutions) {
	  print MUT  $ppos," ",$translation{$codon}," ",$aa,"\n";
	}
      }
    }
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
  unless ($seqfile) {
    print STDERR "$prog error: no sequence file\nAborting.\n";
    exit(6);
  }
  $parameters.="\t--seqfile=$seqfile\n";
  if ($seqformat) {
    $parameters.="\t--seqformat=$seqformat\n";
  }
  unless ($firstaa) {
    $firstaa=1;
  }
  $parameters.="\t--first_aa=$firstaa\n";
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
  if ($exclude_stop) {
    $parameters.="\t--exclude_stop\n";
  }
  if ($prep_file) {
    if ($mutfile || $resfile || $res) {
      die "The flag --prep-mut-file is incompatible with --mutfile, --resfile or --res\n";
    }
    else {
      $parameters.="\t--prep-mut-file = $prep_file\n";
    }
  }
  print $parameters,"\n";
}

sub print_help_message {
    my $help=
"The program $prog reads an mRNA sequence and calculates its mutational \n".
"susbsitution landscape, i.e., all the possible substitutions that involve a \n".
"Given number of changes in the sequence (default: just one change).\n\n".
"For example, if the input sequence is:\n".
">Example\n".
"AUGACGCCC\n".
"The original translated sequence is:\n".
"MTP\n".
"The output peptide seqs will be:\n".
"MTP, LTP, VTP, TTP, KTP, RTP, ITP, MSP, MPP, MAP, MMP, MKP, MRP, MTS, MTT, MTA,\n".
"MTL, MTH, MTR\n";
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
