#!/usr/bin/perl -w

#                       #
# calc_substitutions.pl #
#                       #

# Generate a hash of all possible single genetic code subsitutions

#ccccccccccccccccccccccccccccccccccccccccccccccc#
# Ran Friedman                                  #
# Computaional Chemistry and Biochemistry Group #
# Department of Chemistry and Biomedicinal      #
# Sciences, Linnaeus University, Kalmar, Sweden #
# email: ran.friedman@lnu.se                    #
#ccccccccccccccccccccccccccccccccccccccccccccccc#

use strict;
use List::MoreUtils qw/ uniq /;

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
		   'GGU' => 'G'    # Glycine
		  );

my $codon;
my $i;
my $char;

foreach $codon (sort keys %translation) {
  my @subs;
  for ($i=0;$i<3;$i++) {
    my $newcodon=$codon;
    foreach $char ('A','C','G','U') {
      next if (substr($newcodon,$i,1) eq $char);
      substr($newcodon,$i,1)=$char;
      if ($translation{$newcodon} ne  $translation{$codon}) {
	push(@subs,$translation{$newcodon});
      }
    }
  }
  print "'$codon' => '";
  foreach $char (uniq (sort @subs)) {
    print "$char ";
  }
  print "',\n";
}
